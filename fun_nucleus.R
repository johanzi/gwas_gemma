

#function upload a tab-delimited file and turn into a large list.
#Provides as argument the path containing the file and the type of structure 
#contained in the title of the file ("nuclei" or "CCs")
import.large.list <- function(path, pattern){
  #Create a large_list containing the Summary files for all accessions
  ##Get the name of the files in the directory of interest (files retrieved using bash script)
  names.files <- list.files(path=path, pattern=pattern, full.names=TRUE)
  
  ##Create the large_list using read.delim function to import the files
  list1 <- lapply(names.files, read.delim)
  
  ##Get the names of each files for the dataframe labelling
  names(list1) <- list.files(path=path, pattern=pattern, full.names=FALSE)
  
  #Return the large list generated
  return(list1)
}


# Function returning a dataframe with a new column entitled "Accession", If argument "name" is provided, the date will be removed
# from accession name. If not, the date is included
get.df <- function(large.list, name_mode='default'){
  
  #Tips from http://stackoverflow.com/questions/22689301/how-use-read-table-in-a-for-loop
  
  #Turn the large list into a dataframe (combining by row)
  df <- do.call(rbind.data.frame, large.list)
  
  #Replace blanks by NA values in the dataframe
  df[df==""] <- NA
  
  #Get name of accession for each line from the "Label" column of the dataframe df
  Label <- df$Label
  
  #Convert the list into character to allow the parsing
  Label <- as.character(Label)
  
  #Combine the list by rows
  Label <- rbind(Label)
  
  #Parse the "Label" so that the number of the image + the extension are removed
  # If argument name_mode="name", return only the accession name, otherwise, return name and date
  if(name_mode == "name"){
    Label <- (strsplit(Label, "_"))
  } else {
    Label <- (strsplit(Label, "(-|_)[0-9]{1,}.tif$"))
  }
  
  #Get all first element of the nested lists created (which should contain the accession name)
  accession <- lapply(Label, `[[`, 1)
  
  #add a new column in the dataframe df called "Accession" which contains all 
  #parsed names defined in acc
  df$Accession <- accession
  
  #Put the column "Accession" at the first position (choose last column (Accession)), 
  #New pipeline with nuclei results generated 47 columns
  df <- subset(df, select=c(Accession, 1:47))
  
  #Unlist column "Accession" 
  df$Accession <- unlist(df$Accession)
  
  #Turn the column "Accession" into a factor
  df$Accession <- as.factor(df$Accession)
  
  #Return the dataframe
  return(df)
}

#Summarize nuclei variables for each accession
summarize.nuclei <- function(df){
  
  #Mean the different variables and select only the columns for nuclear characteristics
  
  df_nuclei <- subset(df, select=c(1:26))
  
  df_CC <- subset(df, select=c(1,2, 28:48))
  
  df_mean_nuclei <- aggregate(df_nuclei[3:26], list(df$Accession), mean)
  names(df_mean_nuclei)[1] <- "Accession"
  
  df_mean_CC <- aggregate(df_CC[3:23], list(df$Accession), mean, na.rm = TRUE)
  names(df_mean_CC)[1] <- "Accession"
  
  #Keep accession tag, name, label, and type
  df_type <- subset(df, select=c(49,1,27))
  
  #Get HX values, 1 value per accession
  df_hx <- hx.calculator(df_type)
  
  #Merge dataframes
  df_merged <- merge(df_mean_nuclei, df_mean_CC, by="Accession")
  df_merged <- merge(df_merged, df_hx, by="Accession")
  
  #Add tag
  df_tag <- levels(df$tg_ecotypeid)
  df_accession <- levels(df$Accession)
  df_name <- cbind(df_tag, df_accession)
  df_name <- as.data.frame(df_name)
  names(df_name)[1] <- "tg_ecotypeid"
  names(df_name)[2] <- "Accession"
  
  df_merged <- merge(df_merged, df_name, by="Accession")
  
  #Reorder columns
  df_merged <- subset(df_merged, select=c(1, 48, 2:47))
  
  return(df_merged)
  
}



hx.calculator <- function(df){
  #Hx all accessions together
  general_hx <- sum(df$Type=="compact")/sum(df$Type!="")
  
  #Summarize the value of Type for each accession (creates a table)
  table_type <- (table(df$Accession, df$Type))
  
  #Convert in df
  df_table_type <- as.data.frame.matrix(table_type)
  
  #Convert row names into a column "Accession"
  df_table_type <- tibble::rownames_to_column(df_table_type, "Accession")
  
  #Add a 4th column with the ratio 'compact nuclei/all nuclei'
  df_table_type[4] <- df_table_type[2]/(df_table_type[3]+df_table_type[2])
  
  #Keep only the first and last column
  df_hx <- df_table_type[c(1,4)]
  
  #Rename the 2 columns of the dataframe
  colnames(df_hx) <- c("Accession", "HX")
  
  return(df_hx)
}

#Print a box plot of the chosen dataframe and variable. Title should be added.
ggplot.boxplot <- function(df, group, variable, title){
  ggplot(data = df, aes_string(x=group, y=variable)) + 
    ggtitle(title) + #Title
    theme(plot.title = element_text(hjust = 0.5)) + #Center the title
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_boxplot() #Type of plot (box plot)
}
