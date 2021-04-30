# Function to merge two tables - one with sampling weight and one without, and to write into a csv
# table input is the object created by print
# tbl_nonweight<-print(CreateTableOne())
# tbl_weight<-print(svyCreateTableOne())

merge_table1<-function(tbl_nonweight,tbl_weight,out_path){
  options(stringsAsFactors = FALSE) # make sure strings are not treated as factors in the table
  tbl_weight<-as.data.frame(tbl_weight, stringsAsFactors = FALSE)
  tbl_nonweight<-as.data.frame(tbl_nonweight, stringsAsFactors = FALSE)
  new_tbl<-tbl_weight # use the weighted table as base (including the p, test, missing)
  rownames(new_tbl)<-rownames(tbl_nonweight)
  for (i in 1:(ncol(tbl_weight)-3)){ # loop over each strata
    for (j in 1:nrow(tbl_weight)){ # loop over each row
      if (stringr::str_detect(rownames(tbl_weight)[j],"mean..SD")){ # if "mean  (SD)" is detected in the rowname
        new_tbl[j,i]<-new_tbl[j,i] # then don't modify the weighted table
      } else if (stringr::str_detect(new_tbl[j,i],"^\\s+$")) { # if the cell is empty (only contains white spaces)
        new_tbl[j,i]<-new_tbl[j,i]
      } else { 
        if (str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2]==""){ # if no percentage is calculated (for the "n" row)
          new_tbl[j,i]<-tbl_nonweight[j,i] # then paste only the nonweighted count number
        } else { # otherwise stitch the count from the nonweighted and the percentage from the weighted together
          new_tbl[j,i]<-paste0(str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,1],"(",str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2])
        }
        
      }
    }
  }
  write.csv(new_tbl, file = out_path)
}