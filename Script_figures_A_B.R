# 8/11/2023
# Lucile Massenet-Regad


### Figure A
db_old=readxl::read_excel("~/Downloads/41467_2021_21244_MOESM4_ESM.xlsx", guess_max = 3000) # provided as supplementary material in original article - ICELLNET V1
db_old$`Ligand 3`=NA
db_old$`Ligand 4`=NA
db_old$`Receptor 4`=NA
db_old$`Receptor 5`=NA
db_old$Interaction_name=name.lr.couple(db=db_old, type = "Family")[,1]

db=readxl::read_excel("~/Documents/ICELLNET/Databases/DB_ICELLNET_20230726_classified.xlsx", guess_max = 3000)
table(db$Family, db$Subfamily, useNA = "always")
db$Interaction_name=name.lr.couple(db=db)

db_comp=comparison.lr.couple(db=db) # in V2
db_comp2=comparison.lr.couple(db=db_old) # in V1
length(intersect(db_comp, db_old$Interaction_name)) #in both version
length(intersect(db_comp2, db$Interaction_name)) 

setdiff(db_old$Interaction_name, db_comp)# corrected / excluded in V2 versus V1
db$Update="V2"
db$Update[which(db$Interaction_name %in% intersect(db_comp2, db$Interaction_name))]="V1"

df=as.data.frame(table(db$Family, db$Update, useNA = "ifany"))
head(df)
df$Var1=as.character(df$Var1)
df$Var1[is.na(df$Var1)]="Other"
df$Var1 = factor(df$Var1, levels = df$Var1[order(c(df$Freq[which(df$Var2=="V2")] + df$Freq[which(df$Var2=="V1")]), decreasing = TRUE)])
ggplot(data=df, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity", color="black") + coord_flip() +
  theme_bw() + xlab("Family") + ylab("Number of interactions") +
  scale_fill_manual(values = c("lightgrey", "#0073C2FF"))


### Number for Figure B - Need to run Databases_comparison.R before
Cellphone_DB=db.new
Cellchat_DB=db.new2
head(DB) # ICELLNET database


# ----COMPARISON CellPHONEDB WITH ICELLNET ----
all_names=comparison.lr.couple(db=DB)
length(intersect(all_names, Cellphone_DB$Interaction_name) %>% unique()) # ICELLNET name found in CellPhoneDB 
length(setdiff(Cellphone_DB$Interaction_name,all_names ) )

#symmetry ICELLNET vs CellphoneDB
all_names2=comparison.lr.couple(db=Cellphone_DB)
length(intersect(all_names2, DB$Interaction_name) %>% unique()) # CellphoneDB name included in ICELLNET
length(setdiff(DB$Interaction_name,all_names2))

length(intersect(intersect(all_names, Cellphone_DB$Interaction_name), intersect(all_names2, DB$Interaction_name))) #900 symmetrical names found

sort(setdiff( intersect(all_names2, DB$Interaction_name), intersect(all_names, Cellphone_DB$Interaction_name))) # From CellPhoneDB, found only in ICELLNET list
sort(setdiff(intersect(all_names, Cellphone_DB$Interaction_name), intersect(all_names2, DB$Interaction_name))) #from ICELLNET , found only in CellphoneDB list -- 1 more found here -> IFNL3 / IFNLR1 + IL10RB because included 2 times in CellPhoneDB data

all(setdiff(intersect(all_names, Cellphone_DB$Interaction_name), intersect(all_names2, DB$Interaction_name)) %in% all_names2) #from ICELLNET , found only in CellphoneDB list
all(setdiff( intersect(all_names2, DB$Interaction_name), intersect(all_names, Cellphone_DB$Interaction_name))  %in% all_names) # From CellPhoneDB, found only in ICELLNET list 
###CONCLUSION : do both symmetric operations and take the smallest number + check that list present in the combination derived from other interactions lists --> means there are doublons


# ----COMPARISON CELLCHAT WITH ICELLNET ----
all_names=comparison.lr.couple(db=DB) 
length(intersect(all_names, Cellchat_DB$Interaction_name)) 
length(setdiff(Cellchat_DB$Interaction_name, all_names)) 

all_names2=comparison.lr.couple(db=Cellchat_DB) 
length(intersect(all_names2, DB$Interaction_name)  %>% unique()) 
length(setdiff( DB$Interaction_name, all_names2 )) 

#Same story than for previous DB -- comparison
length(intersect(intersect(all_names, Cellchat_DB$Interaction_name), intersect(all_names2, DB$Interaction_name))) #719 symmetrical names found

sort(setdiff( intersect(all_names2, DB$Interaction_name), intersect(all_names, Cellchat_DB$Interaction_name))) # From CellChat, found only in ICELLNET list
sort(setdiff(intersect(all_names, Cellchat_DB$Interaction_name), intersect(all_names2, DB$Interaction_name))) #from ICELLNET , found only in CellChat list 

all(setdiff(intersect(all_names, Cellchat_DB$Interaction_name), intersect(all_names2, DB$Interaction_name)) %in% all_names2) 
all(setdiff( intersect(all_names2, DB$Interaction_name), intersect(all_names, Cellchat_DB$Interaction_name))  %in% all_names) 
###CONCLUSION : do symmetry and take the smallest number + check that list present in the combination derived from other interactions lists --> means there are doublons


# ----COMPARISON CELLCHAT WITH CELLPHONEDB ---- 
all_names_CellPhone=comparison.lr.couple(db=Cellphone_DB) 
length(intersect(all_names_CellPhone, Cellchat_DB$Interaction_name) %>% unique()) #Number of interaction found in CellPhoneDB
length(setdiff(Cellchat_DB$Interaction_name, all_names_CellPhone) %>% unique())

all_names_CellChat=comparison.lr.couple(db=Cellchat_DB) 
length(intersect(all_names_CellChat, Cellphone_DB$Interaction_name)  %>% unique() ) # Check that doublons are the same (8 similar + 1 another -> CD96/NECTIN1)
length(setdiff(Cellphone_DB$Interaction_name, all_names_CellChat) %>% unique())

all(setdiff(intersect(all_names_CellPhone, Cellchat_DB$Interaction_name), intersect(all_names_CellChat, Cellphone_DB$Interaction_name)) %in% all_names_CellChat) 
all(setdiff( intersect(all_names_CellChat, Cellphone_DB$Interaction_name), intersect(all_names_CellPhone, Cellchat_DB$Interaction_name))  %in% all_names_CellPhone) 

setdiff(intersect(all_names_CellPhone, Cellchat_DB$Interaction_name), intersect(all_names_CellChat, Cellphone_DB$Interaction_name)) %>% sort() %>% as.data.frame() %>% dim()
setdiff( intersect(all_names_CellChat, Cellphone_DB$Interaction_name), intersect(all_names_CellPhone, Cellchat_DB$Interaction_name))  %>% sort() %>% as.data.frame() %>% dim()



# ---- IN THE 3 DATABASES # - -----

length(intersect(all_names, intersect(all_names_CellChat, Cellphone_DB$Interaction_name))) #all names that can be created from ICELLNET interactions 
length(intersect(all_names_CellChat, intersect(all_names_CellPhone, DB$Interaction_name)))
length(intersect(all_names_CellPhone, intersect(all_names_CellChat, DB$Interaction_name)))
length(intersect(all_names_CellPhone, intersect(all_names, Cellchat_DB$Interaction_name)))
length(intersect(all_names, intersect(all_names_CellPhone, Cellchat_DB$Interaction_name)))  
in_three=intersect(all_names, intersect(all_names_CellPhone, Cellchat_DB$Interaction_name))

intersect(intersect(all_names, intersect(all_names_CellPhone, Cellchat_DB$Interaction_name)) ,  intersect(all_names, intersect(all_names_CellChat, Cellphone_DB$Interaction_name)))

all( setdiff(intersect(all_names, intersect(all_names_CellPhone, Cellchat_DB$Interaction_name)) ,  intersect(all_names, intersect(all_names_CellChat, Cellphone_DB$Interaction_name))) %in% all_names_CellChat)
all( setdiff(intersect(all_names, intersect(all_names_CellChat, Cellphone_DB$Interaction_name)), intersect(all_names, intersect(all_names_CellPhone, Cellchat_DB$Interaction_name))) %in% all_names_CellPhone)

# IN CELLCHAT AND CELLPHONEDB but not in ICELLNET 
all(setdiff(intersect(all_names_CellPhone, Cellchat_DB$Interaction_name), intersect(all_names_CellChat, Cellphone_DB$Interaction_name)) %in% all_names_CellChat) 

setdiff(setdiff(intersect(all_names_CellChat, Cellphone_DB$Interaction_name), in_three), all_names) ##343 interactions not included at all in ICELLNET but in the two others DB



