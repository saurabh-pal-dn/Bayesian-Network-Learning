# Bayesian-Network-Learning

Inputs:
alarm.bif(Bayesian Inference file) which has relation of components
gold_alarm.bif which ahs the correct probabilities
records.txt file that has records of 11000 patients but each records missing one entry
BayesNet.png which shows the newtork in pictoral form

Output:
solved_alarm.bif which matches golden_alarm.bif

Algorithm used:
Expectation-Maximization

As we don't have proper information we first expect the probability and then try to maximize our expectation
I used laplace smoothing
Implemeted an epsilon of 0.0005 for computing convergence.

command:
python main.py
