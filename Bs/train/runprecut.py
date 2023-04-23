#! /usr/bin/python3

from subprocess import run, PIPE

inputs = "/home/tasheng/braa/Unskimmed/NewOfficialMC/BsMC.root"
inputb = "/home/tasheng/braa/Unskimmed/BsData.root"
output = "BsMCCut.root"

rootarg = f'precut.cc+("{inputs}", "{inputb}", "1", "1", "{output}")'
# print(rootarg)
run(['root', '-q', '-b', '-l', rootarg ])
