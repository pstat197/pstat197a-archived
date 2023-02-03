library(reticulate)
path_to_python <- 'c:/users/trdea/minconda3/python.exe'
virtualenv_create("r-reticulate", python = path_to_python)

library(keras)
install_keras()

library(tensorflow)
tf$constant("Hello Tensorflow!")
