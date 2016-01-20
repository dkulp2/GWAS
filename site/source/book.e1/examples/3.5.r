# Example 3.5 (Population substructure and LD):
ObsCount <- matrix(c(136,64,64,136),2)
ObsCount
ExpCount <- chisq.test(ObsCount)$expected
ExpCount