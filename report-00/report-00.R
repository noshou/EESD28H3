library(datasets)
library(tidyverse)

# load iris dataset
data(iris)

# Q1 #

# get column names
col_names <- colnames(iris)
print(col_names)

# map each column data type (class) to a list
data_types <- iris %>% map(class)
print(data_types)

# Q2 #
print(nrow(iris))
print(ncol(iris))

# Q3 #
width <- iris$Sepal.Width

# Q4 #
print(width[100])

# Q5 #
print(width[length(width)])

# Q6 #
iris_rows_10_to_20 <- iris %>% slice(10:20)
print(iris_rows_10_to_20)

# Q7 #

# method 1; just select rows
method_1 <- iris_rows_10_to_20 %>% select(Species, Petal.Width, Petal.Length)
print(method_1)

# method 2; "pull" out vectors and concatenate together
method_2 <- data.frame(
	Species      =	iris_rows_10_to_20 %>% pull(Species), 
	Petal.Width  = 	iris_rows_10_to_20 %>% pull(Petal.Width),
	Petal.Length =  iris_rows_10_to_20 %>% pull(Petal.Length)
)
print(method_2)

# Q8 #
iris_row_slice <- iris %>% slice(1:10,20,100)
print(iris_row_slice)

# Q9 #
method_1 <- iris_row_slice$Sepal.Length[1]
print(method_1)
method_2 <- (iris_row_slice %>% pull(Sepal.Length))[1]
print(method_2)
method_3 <- first(iris_row_slice$Sepal.Length)
print(method_3)

# Q10 #
# option a would work, since we concatenate 1,2,3 into a vector [1,2,3]
# which represents row 1,2,3 of the Sepal.Length column.
# option b would not work since 
