using SLHQuantumSystems

SQZ = squeezing_cavity("A")
FCV = cavity("B")

SYS1 = concatenate([SQZ,FCV],"sys")
SYS2 = feedbackreduce(SYS1,"A_out","B_in")

SS = StateSpace(SYS2)
