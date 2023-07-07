
import os
import uproot
import optuna
from   optuna.trial import TrialState
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data
from   torchvision import datasets
from   torchvision import transforms
import numpy as np
import random


DEVICE = torch.device("cpu")
BATCHSIZE = 128
CLASSES = 10
DIR = os.getcwd()
EPOCHS = 10
N_TRAIN_EXAMPLES = BATCHSIZE * 30
N_VALID_EXAMPLES = BATCHSIZE * 10


def define_model(trial):
    # We optimize the number of layers, hidden units and dropout ratio in each layer.
    n_layers = trial.suggest_int("n_layers", 1, 3)
    layers = []

    in_features = 28 * 28
    for i in range(n_layers):
        out_features = trial.suggest_int("n_units_l{}".format(i), 4, 128)
        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = trial.suggest_float("dropout_l{}".format(i), 0.2, 0.5)
        layers.append(nn.Dropout(p))

        in_features = out_features
    layers.append(nn.Linear(in_features, CLASSES))
    layers.append(nn.LogSoftmax(dim=1))

    return nn.Sequential(*layers)


def get_root(stages,ptmin,ptmax):
    # Load FashionMNIST dataset.
    fileS = uproot.open("/lstore/cms/simao/sample/BPMC_3_60.root")
    fileB = uproot.open("/lstore/cms/simao/sample/BPData_3_60.root")

    treeS = fileS["Bfinder/ntKp"]
    treeShl = fileS["hltanalysis/HltTree"]
    treeShi = fileS["hiEvtAnalyzer/HiTree"]
    treeSskim = fileS["skimanalysis/HltTree"]
    treeB = fileB["Bfinder/ntKp"]
    treeBhl = fileB["hltanalysis/HltTree"]
    treeBhi = fileB["hiEvtAnalyzer/HiTree"]
    treeBskim = fileB["skimanalysis/HltTree"]

    treeS_full=np.concatenate((treeS,treeShl,treeShi,treeSskim),1)
    treeB_full=np.concatenate((treeB,treeBhl,treeBhi,treeBskim),1)

    cut="((pPAprimaryVertexFilter == 1) & (pBeamScrapingFilter == 1) & (HLT_HIL1DoubleMu0_v1 == 1))  &  ((Bmu1isTriggered == 1) & (Bmu2isTriggered == 1) ) & ((Btrk1Pt > 0.5 && Bchi2cl > 0.05) & (BsvpvDistance/BsvpvDisErr > 2.0) & (Bpt > 2) & (abs(Btrk1Eta-0.0) < 2.4)  & ((TMath::Abs(By)<2.4) & (TMath::Abs(Bmumumass-3.096916)<0.15) & (((abs(Bmu1eta)<1.2) & (Bmu1pt>3.5)) || ((abs(Bmu1eta)>1.2) & (abs(Bmu1eta)<2.1) & Bmu1pt>(5.47-1.89*abs(Bmu1eta))) || (abs(Bmu1eta)>2.1 & abs(Bmu1eta)<2.4 & Bmu1pt>1.5)) & ((abs(Bmu2eta)<1.2 & Bmu2pt>3.5)||(abs(Bmu2eta)>1.2 & abs(Bmu2eta)<2.1 & Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1 & abs(Bmu2eta)<2.4 & Bmu2pt>1.5)) & Bmu1TMOneStationTight & Bmu2TMOneStationTight & Bmu1InPixelLayer>0 & (Bmu1InPixelLayer+Bmu1InStripLayer)>5 & Bmu2InPixelLayer>0 & (Bmu2InPixelLayer+Bmu2InStripLayer)>5 & Bmu1dxyPV<0.3 & Bmu2dxyPV<0.3 & Bmu1dzPV<20 & Bmu2dzPV<20 & Bmu1isTrackerMuon & Bmu2isTrackerMuon & Bmu1isGlobalMuon & Bmu2isGlobalMuon & Btrk1highPurity & abs(Btrk1Eta)<2.4 & Btrk1Pt>0.5)   &  (Btrk1PixelHit + Btrk1StripHit > 10)  &   (Btrk1PtErr/Btrk1Pt < 0.1) &  Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18    &  (abs(PVz)<15))"
    cuts="%s & Bgen==23333 & Bpt>%f & Bpt<%f " % (cut, ptmin, ptmax)
    cutb="%s & ((Bmass - 5.27929 ) > 0.25 &  (Bmass - 5.27929) < 0.30) & Bpt>%f & Bpt<%f" % (cut, ptmin, ptmax)
    
    signal=treeS_full.arrays(cut=cuts,aliases={"Trk1DCAz":"abs(Btrk1Dz1/Btrk1DzError1)","Trk2DCAz":"abs(Btrk2Dz1/Btrk2DzError1)","Trk1DCAxy":"abs(Btrk1Dxy1/Btrk1DxyError1)","Trk2DCAxy":"abs(Btrk2Dxy1/Btrk2DxyError1)","MassDis":"abs(Btktkmass-1.019455)","dls":"BsvpvDistance/BsvpvDisErr","dls2D":"Bd0"})
    background=treeB_full.arrays(cut=cutb,aliases={"Trk1DCAz":"abs(Btrk1Dz1/Btrk1DzError1)","Trk2DCAz":"abs(Btrk2Dz1/Btrk2DzError1)","Trk1DCAxy":"abs(Btrk1Dxy1/Btrk1DxyError1)","Trk2DCAxy":"abs(Btrk2Dxy1/Btrk2DxyError1)","MassDis":"abs(Btktkmass-1.019455)","dls":"BsvpvDistance/BsvpvDisErr","dls2D":"Bd0"})

    nsignal=len(signal["Bmass"])
    nback=len(background["Bmass"])
    nevents=len(nsignal+nback)
    x=np,zeros([nevents,len(stages)])
    y=np.zeros(nevents)
    y[:nsignal]=1
    for i,j in enumerate(stages):
        x[:nsignal,:]=signal[j]
        x[nsignal:,:]=background[j]

    

    data= {"train": (train_X, data["ytrain"]),
            "test": (test_X, data["ytest"])}



    train_loader = torch.utils.data.DataLoader(
        datasets.FashionMNIST(DIR, train=True, download=True, transform=transforms.ToTensor()),
        batch_size=BATCHSIZE,
        shuffle=True,
    )
    valid_loader = torch.utils.data.DataLoader(
        datasets.FashionMNIST(DIR, train=False, transform=transforms.ToTensor()),
        batch_size=BATCHSIZE,
        shuffle=True,
    )

    return train_loader, valid_loader


def objective(trial):
    # Generate the model.
    model = define_model(trial).to(DEVICE)

    # Generate the optimizers.
    optimizer_name = trial.suggest_categorical("optimizer", ["Adam", "RMSprop", "SGD"])
    lr = trial.suggest_float("lr", 1e-5, 1e-1, log=True)
    optimizer = getattr(optim, optimizer_name)(model.parameters(), lr=lr)

    # Get the FashionMNIST dataset.
    train_loader, valid_loader = get_root()

    # Training of the model.
    for epoch in range(EPOCHS):
        model.train()
        for batch_idx, (data, target) in enumerate(train_loader):
            # Limiting training data for faster epochs.
            if batch_idx * BATCHSIZE >= N_TRAIN_EXAMPLES:
                break

            data, target = data.view(data.size(0), -1).to(DEVICE), target.to(DEVICE)

            optimizer.zero_grad()
            output = model(data)
            loss = F.nll_loss(output, target)
            loss.backward()
            optimizer.step()

        # Validation of the model.
        model.eval()
        correct = 0
        with torch.no_grad():
            for batch_idx, (data, target) in enumerate(valid_loader):
                # Limiting validation data.
                if batch_idx * BATCHSIZE >= N_VALID_EXAMPLES:
                    break
                data, target = data.view(data.size(0), -1).to(DEVICE), target.to(DEVICE)
                output = model(data)
                # Get the index of the max log-probability.
                pred = output.argmax(dim=1, keepdim=True)
                correct += pred.eq(target.view_as(pred)).sum().item()

        accuracy = correct / min(len(valid_loader.dataset), N_VALID_EXAMPLES)

        trial.report(accuracy, epoch)

        # Handle pruning based on the intermediate value.
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return accuracy


if __name__ == "__main__":
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=100, timeout=600)

    pruned_trials = study.get_trials(deepcopy=False, states=[TrialState.PRUNED])
    complete_trials = study.get_trials(deepcopy=False, states=[TrialState.COMPLETE])

    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    print("Best trial:")
    trial = study.best_trial

    print("  Value: ", trial.value)

    print("  Params: ")
    for key, value in trial.params.items():
        print("    {}: {}".format(key, value))