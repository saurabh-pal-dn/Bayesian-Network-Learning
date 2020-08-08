# Now we will use Expectation Maximiztion to read the stuff from records.txt and as they are in large number we will basically have a normal 
# distributon for them. 
# Once we do this we will do some magic(Figure out)
# Then its simply using these values to do MCMC Gibbs sampling for the network to get the probabilities and then feed them into the CPT
import time
import runner as rn
import BayesianNetwork 

step0=time.time()
print("Initializing.....")
alarm,df,missingIndex=rn.setupNetwork("Alarm.bif","records.txt")
step1=time.time()
print("Initializing time: (%ss)" % round(step1-step0,5))
# print("Printing Old network")
# BayesianNetwork.printNetwork(alarm)
print("Learning parameters.....")
step2=time.time()
newalarm=rn.ExpMax(df,alarm,missingIndex)
print("Learning parameters time: (%ss)" % round(step2-step1,5))
print("Printing new network")
BayesianNetwork.printNetwork(newalarm)