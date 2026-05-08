#Read in test data 
Test_data_set <- read.csv("Test data set.csv")
head(Test_data_set)

# Make a table of GSE 
table(Test_data_set$GSE)
#colouring bars based on GSE
Test_data_set$GSE <- as.factor(Test_data_set$GSE)
library(ggplot2)

ggplot(Test_data_set, aes(x = Age, fill = GSE)) +
  geom_histogram(binwidth = 1, position = "stack") +
  scale_fill_manual(values = c("GSE185920" = "skyblue", 
                               "GSE185445" = "pink", 
                               "GSE149318" = "green")) +
  ylim(0, 100) +
  labs(title = 'Distribution of ages in test data set', 
       x = "Age", y = "Frequency", fill = "GSE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
