# Reproducible research: version control and R

## For questions 1, 2 and 3: 

## 4. Simulating a random walk
```{r}
#install.packages("ggplot2")
#install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}
```
A dataframe is created containing 3 variables: x, y and time. The code simulates a 'random walk' according to x- and y-coordinates. A for loop is used and creates an iteration, whereby the size of the step and the angle is specified. The subsequent step (in terms of the angle/direction) is determined from the previous step.

```{r}
data1 <- random_walk(500)

plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")
```

```{r}

data2 <- random_walk(500)

plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

grid.arrange(plot1, plot2, ncol=2)
```

Here, we are creating 2 different 'random walk' plots using the ggplot2 package, each with 500 steps, and putting them next to one another for comparison. The for-loop used enables an adjustment to the angle and step size with every step taken. Since this is generated randomly, the graphs are completely different and are predicted to change every time the code is re-run. This is because of the fact that the direction taken with respect to time is generated at random (from the coordinates of the previous step). 
Time is measured here in number of steps, and demonstrated by use of a colour key - the darker the shade of blue the more time has elapsed (and thus the more steps taken).


## 4b. Random seeds

A random seed is a method used in R for generating a pseudorandom number. The number is an integer vector generated with an algorithm, but requires a 'seed' to initialise. Hence, the number produced is pseudorandom because if you know both the seed and the generator, you can predict and reproduce the outcome.

The algorithm random number generator (RNG) mimics the properties of the independent generation of numbers within a distribution in the interval (0,1). Therefore, the random seed is the first number used for this generation of numbers.

Another reason the random seed is useful is because it ensures reproducibility of results, meaning the output of random numbers will be the same in each run.
(information was sourced from r-coder.com)


## 4d. Comparison of the random_walk code with the new brownian_motion_code

<img width="1387" alt="Random walk vs brownian motion" src="https://github.com/lauren-homan/reproducible-research_homework/assets/150149060/e12c3d43-d34d-450f-a050-417b7790e67b">


## 5. How many rows and columns does the dsDNA dataset contain?

There are 13 columns and 33 rows.

## 5b. What transformation can you use to fit a linear model to the data?

 **$`V = \beta L^{\alpha}`$**
```{r}
install.packages("janitor")
install.packages('dplyr')
#I am installing these packages - janitor will be necessary in order to clean the code, and dplyr will help organise the code

library(janitor)
library(dplyr)


viral_data<- read.csv("question-5-data/Cui_etal2014.csv")


viral_clean <- viral_data %>%
  clean_names()
names(viral_clean)
```

I have firstly cleaned the data so that the column headings contain no capitals and are in snake case  

For the transformation, we must log-transform both V and L such that **$`log(V) = {\alpha}log(L)+\beta`$**

This will then enable us to apply a linear model, as demonstrated in the code below:

```
viral_clean <- viral_clean %>%
  mutate(log_virion_volume = log(virion_volume_nm_nm_nm))

viral_clean <- viral_clean %>% 
  mutate(log_genome_length = log(genome_length_kb))


linear_model<- lm (log_virion_volume ~ log_genome_length, data = viral_clean)
summary(linear_model)
```
## 5c. Finding exponent and scaling factor

<img width="455" alt="Screenshot 2023-12-06 at 18 04 59" src="https://github.com/lauren-homan/reproducible-research_homework/assets/150149060/9ee6cd5c-a4d2-42f4-b4ef-1da27e84e132">

As shown in the screenshot:

**$`{\alpha}`$** = 1.515

p-value = 2.28e-10

**$`\beta`$** = e^7.0748 = 1181.81

p-value = 6.44e-10


Due to very small p-values, we can conclude that the values for the exponent and scaling factor within our model and obtained from our data are significant.

#### In the paper:

Exponent (**$`{\alpha}`$**) = 1.43

Scaling factor (**$`\beta`$**) = 2057


The values for the exponent, alpha, are quite closely aligned, however the scaling factor found in the paper is much larger than that of our data. The paper includes confidence intervals for the scaling factor (1185-3571), meaning our value for beta is significantly different to the paper's. Perhaps this is due to differences in sample sizes, which would alter the values derived from the model.


## 5d. The code to make the plot is as shown below

```{r}

ggplot(viral_clean, aes(x = log_genome_length, y = log_virion_volume))+
  geom_point() +
  geom_smooth(method='lm', se = TRUE) +
  theme_bw() +
  labs( x = 'log[Genome length (kb)]',
        y = 'log[Virion volume (nm3)]')

```
## 5e. Estimated volume of dsDNA virus

**$`V = \beta L^{\alpha}`$**

When:

L = 300kb


V = 6698076


## Bonus question

Reproducability and replicability are similar terms used in science, but with slightly different meanings in reference to the scientific method. Reproducability is defined as the ability to obtain the same results when the same experiment is performed by different scientists. Replicability, on the other hand, is when a similar (but not identical) experiment is replicated under different conditions, using different data or slightly modifying the methods.

Therefore, while both processes hinge upon reliability of scientific method, replicability utilises the same methodology with different conditions whereas reproducibility entirely repeats an investigation.


The Git and GitHub platforms were designed to enable sharing, collaboration and integration of data. They enable collaborators to build upon previous work remotely, or indeed, to use other people's work for their own analyses. Collaboration is a central premise within science therefore that these platforms have made the process easier is a huge benefit to scientific progression.

These platforms enable reproducibility and replicability by providing a comprehensive step-by-step guide for other researchers to work on, meaning projects can be managed from different parts of the world if necessary, without losing any significant progress. The ability to observe any edits or changes made to the code is also highly useful, as it is clear what has been changed and why.


With that being said, GitHub is a difficult platform to navigate when using it for the first time. Some of the steps required when editing a repository are not very clear, meaning the user may potentially end up learning how to use the platform but without really understanding what it is they're doing. This is of no benefit, considering its entire purpose is to enable collaboration and data sharing. It is a hugely useful platform, but in order to become used on a large scale it may benefit from detailed tutorials and instructions.

## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   - A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
   - Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
   - Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
   - Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \beta L^{\alpha}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   - Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
   - What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
   - Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
   - Write the code to reproduce the figure shown below. (10 points)

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  - What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
