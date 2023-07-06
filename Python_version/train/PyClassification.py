

import argparse
import torch

import uproot
import torch.nn as nn
from matplotlib import pyplot as plt



class FeedforwardNetwork(nn.Module):
    def __init__(
            self, n_classes, n_features, hidden_size, layers,
            activation_type, dropout, **kwargs):
        
        super().__init__()

        self.first_layer = nn.Linear(n_features, hidden_size)
        self.hidden_layers = nn.ModuleList([nn.Linear(hidden_size, hidden_size) for i in range(layers-1)]) #creates a list which holds modules, letting us automatically create different nn.Linears for the hidden layers.

        self.output_layer = nn.Linear(hidden_size, n_classes)

        if activation_type == "relu": self.activation = nn.ReLU()
        else: self.activation = nn.Tanh()
        drop=nn.Dropout(p=dropout)
        self.order = nn.Sequential(self.first_layer,self.activation,drop)   
        for k in range(layers-1):                                           
            self.order.append(self.hidden_layers[k])
            self.order.append(self.activation)
            self.order.append(drop)
        self.order.append(self.output_layer)                                
        
        
        
    def forward(self, x, **kwargs):
        
        x = self.order(x)                             
        return x


def train_batch(X, y, model, optimizer, criterion, **kwargs):

 
    optimizer.zero_grad()    # reset the grad value
    y2 = model(X)            # predicted scores
    loss = criterion(y2,y)   # needs to be (pred_scores,gold labels), y2 has an extra dimension due to tracking of grad
    loss.backward()          # calculates the grads
    optimizer.step()         # updates weight 
    return loss.item()       # only want the value, giving everything occupies too much memory
    


def predict(model, X):
    
    scores = model(X)  # (n_examples x n_classes)
    predicted_labels = scores.argmax(dim=-1)  # (n_examples)
    return predicted_labels


def evaluate(model, X, y):
   
    model.eval()
    y_hat = predict(model, X)
    n_correct = (y == y_hat).sum().item()
    n_possible = float(y.shape[0])
    model.train()
    return n_correct / n_possible


def plot(epochs, plottable, ylabel='', name=''):
    plt.clf()
    plt.xlabel('Epoch')
    plt.ylabel(ylabel)
    plt.plot(epochs, plottable)
    plt.savefig('%s.pdf' % (name), bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ptmin', type=float)
    parser.add_argument('ptmax', type=float)
    parser.add_argument('-epochs', default=30, type=int)
    parser.add_argument('-batch_size', default=256, type=int)
    parser.add_argument('-hidden_size', type=int, default=128)
    parser.add_argument('-layers', type=int, default=20)
    parser.add_argument('-dropout', type=float, default=0.0)
    parser.add_argument('-stages', type=int, default=[0,2,4,7,8,11])
    parser.add_argument('-activation',
                        choices=['tanh', 'relu'], default='tanh')
    parser.add_argument('-optimizer',
                        choices=['sgd', 'adam'], default='adam')
    opt = parser.parse_args()

    
    data = utils.load_classification_data()
    dataset = utils.ClassificationDataset(data)
    train_dataloader = DataLoader(
        dataset, batch_size=opt.batch_size, shuffle=True)

    dev_X, dev_y = dataset.dev_X, dataset.dev_y
    test_X, test_y = dataset.test_X, dataset.test_y

    n_classes = torch.unique(dataset.y).shape[0]  # 10
    n_feats = dataset.X.shape[1]

    # initialize the model

    model = FeedforwardNetwork(
        n_classes,
        n_feats,
        opt.hidden_size,
        opt.layers,
        opt.activation,
        opt.dropout
    )

    # get an optimizer
    optims = {"adam": torch.optim.Adam, "sgd": torch.optim.SGD}

    optim_cls = optims[opt.optimizer]
    optimizer = optim_cls(
        model.parameters(),
        lr=opt.learning_rate,
        weight_decay=opt.l2_decay)

    # get a loss criterion
    criterion = nn.CrossEntropyLoss()

    # training loop
    epochs = torch.arange(1, opt.epochs + 1)
    train_mean_losses = []
    valid_accs = []
    train_losses = []
    for ii in epochs:
        print('Training epoch {}'.format(ii))
        for X_batch, y_batch in train_dataloader:
            loss = train_batch(
                X_batch, y_batch, model, optimizer, criterion)
            train_losses.append(loss)

        mean_loss = torch.tensor(train_losses).mean().item()
        print('Training loss: %.4f' % (mean_loss))

        train_mean_losses.append(mean_loss)
        valid_accs.append(evaluate(model, dev_X, dev_y))
        print('Valid acc: %.4f' % (valid_accs[-1]))

    print('Final Test acc: %.4f' % (evaluate(model, test_X, test_y)))
    # plot
    if opt.model == "logistic_regression":
        config = "{}-{}".format(opt.learning_rate, opt.optimizer)
    else:
        config = "{}-{}-{}-{}-{}-{}-{}".format(opt.learning_rate, opt.hidden_size, opt.layers, opt.dropout, opt.activation, opt.optimizer, opt.batch_size)

    plot(epochs, train_mean_losses, ylabel='Loss', name='{}-training-loss-{}'.format(opt.model, config))
    plot(epochs, valid_accs, ylabel='Accuracy', name='{}-validation-accuracy-{}'.format(opt.model, config))


if __name__ == '__main__':
    main()
