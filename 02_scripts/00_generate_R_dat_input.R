# addicts dataset
# These data comprise the times in
# days spent by heroin addicts from entry to departure
# from one of two methadone clinics. There are two further
# covariates, namely, prison record and methadone
# dose, believed to affect the survival times. 
# Column 1: Subject ID
# Column 2: Clinic (1 or 2)
# Column 3: Survival status (0 = censored, 1 = departed from clinic)
# Column 4: Survival time in days
# Column 5: Prison record (0 ¼ none, 1 ¼ any)
# Column 6: Methadone dose (mg/day)

addicts <- read.table(file = "01_data/addicts.dat", skip = 19)
names(addicts) <- c('id', 'clinic', 'status', 'survt', 'prison', 'dose')
saveRDS(addicts, file = "01_data/addicts.rds")

# anderson dataset
# consists of remission survival times on 42 leukemia patients,
# half of whom get a certain new treatment therapy and
# the other half of whom get a standard treatment
# therapy. The exposure variable of interest is treatment
# status (Rx ¼ 0 if new treatment, Rx ¼ 1 if standard
#         treatment). Two other variables for control as potential
# confounders are log white blood cell count (i.e., logwbc) and sex. 
# Failure status is defined by the relapse
# variable (0 if censored, 1 if failure).
anderson <- read.table(file = "01_data/anderson.dat", as.is = TRUE, header = TRUE)
colnames(anderson) <- c("Survt", "Relapse", "Sex", "logWBC", "Rx")
anderson$id <- rownames(anderson)
anderson <- anderson[,c("id", "Survt", "Relapse", "Sex", "logWBC", "Rx")]
saveRDS(anderson, file = "01_data/anderson.rds")

# bladder
# The bladder cancer dataset contains recurrent event outcome
# information for eighty-six cancer patients followed
# for the recurrence of bladder cancer tumor after transurethral
# surgical excision (Byar and Green 1980). The exposure
# of interest is the effect of the drug treatment of thiotepa.
# Control variables are the initial number and initial size of
# tumors. The data layout is suitable for a counting processes
# approach. The variables are defined as follows:
# ID – Patient ID (may have multiple observations for the same subject)
# EVENT – Indicates whether the patient had a tumor (coded 1) or not (coded 0)
# INTERVAL – A counting number representing the order of the time interval for a given subject (coded 1 for the
#                                        subject’s first time interval, coded 2 for a subject’s
#                                        second time interval, etc.)
# START – The starting time (in months) for each interval
# STOP – The time of event (in months) or censorship for each interval
# TX – Treatment status (coded 1 for treatment with thiotepa and 0 for the placebo)
# NUM – The initial number of tumors
# SIZE – The initial size (in centimeters) of the tumor

bladder <- read.table(file = "01_data/bladder.dat", skip = 2)
colnames(bladder) <- c( "id", "event", "interval", "inttime",
                     "start", "stop", "tx", "num", "size")
saveRDS(bladder, file = "01_data/bladder.rds")

# vets
# Survival data for 137 patients from Veteran's Administration Lung Cancer Trial.
# Data from Kalbfleisch, J., and Prentice, R., The Statistical Analysis of 
# Failure Time Data, John Wiley and Sons, New York, 1980.
# 
# Column 1 = treatment (1 = standard, 2 = test)
# Column 2 = cell type 1 (1 = large, 0 = other)
# Column 3 = cell type 2 (1 = adeno, 0 = other)
# Column 4 = cell type 3 (1 = small, 0 = other)
# Column 5 = cell type 4 (1 = squamous, 0 = other)
# Column 6 = survival time (days)
# Column 7 = performance status (0 = worst, ..., 100 = best)
# Column 8 = disease duration (months)
# Column 9 = age (years)
# Column 10 = prior therapy (0 = none, 10 = some)
# Column 11 = status (0= censored, 1 = Died)

library(foreign)
vets <- read.dta("http://web1.sph.emory.edu/dkleinb/allDatasets/surv2datasets/vets.dta")
saveRDS(vets, file = "01_data/vets.rds")
