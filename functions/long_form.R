# function to convert wide-form data to long-form
long_form <- function(dataTable,
                      data.cols,
                      category.cols,
                      colnames = c("category", "data")){
  
 rep.cols<-dim(category.cols)[2]
 
 long_data<-cbind(data.cols[rep(rownames(data.cols), each=rep.cols),],
                  category = rep(colnames(category.cols), dim(category.cols)[1]),
                  data = c(t(category.cols)))
 colnames(long_data) = c(data.cols, colnames)
 long_data <- as.data.frame(long_data)
 long_data$data <- as.numeric(as.character(long_data$data))
 
  # remove zero abundances (they can be re-instated via table() if
  # necessary).
  long_data <- long_data[long_data$data > 0, ]
  
  return(long_data)
}