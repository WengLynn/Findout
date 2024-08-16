findout_APO =function (filename, key_pattern) 
{
  library(tidyr)
  library(stringr)
  ref_file = "reffa.fasta"
  cod = paste0("grep -A1 '", key_pattern, "' ", filename, " > ", 
               ref_file)
  system(cod)
  reffa = read.table(ref_file, sep = "\n")
  reffa = as.character(reffa[2, 1])
  mmref = unlist(strsplit(reffa, ""))
  mmcount = length(mmref)
  coln = c(1:mmcount)
  write.table(t(as.data.frame(c("Type", "loc", "NNfrom", "NNto","type_count","all_snp_count","freq","combined_freq","seqname"))), paste0(filename,"_vs_",key_pattern,"_APO_info.txt"), col.names = F, row.names = F, 
              sep = "\t", quote = F, append = F)
  lineCnt = 0
  con <- file(filename, "r")
  while (1) {
    oneline = readLines(con, n = 1)
    if (length(oneline) == 0) {
      break
    }
    if (length(grep(">", oneline)) == 1) {
      name = oneline
      namestr = unlist(str_split(name, "\\|"))
    }
    if (length(grep(">", oneline)) == 0) {
      mm = unlist(strsplit(oneline, ""))
      mmcount = ifelse(mm == "N", 1, 0)
      mmcount2 = ifelse(mm %in% c("-", "A", "T", "C", "G", "N"), 0, 1)
      mmsum = sum(mmcount) #count N
      mmsum2 = sum(mmcount2) #count specific characteristic
      #  if (mmsum < 15 & mmsum2 < 50) {
      #   write.table(lineCnt, "QC_line.txt", sep = "\n", quote = F, col.names = F, row.names = F, na = "", append = T)
      lineCnt = lineCnt + 1
      if (lineCnt%%10000 == 0) {
        print(lineCnt)
      }
      
      mmc = mm
      mmref_temp = mmref
      
      if (length(mm) != length(mmref)) {
        mmref_temp = mmref[1:length(mm)]
      }
      mmref_temp[!mm %in% c("-", "A", "T", "C", "G")] = mm[!mm %in% c("-", "A", "T", "C", "G")]
      
      loc_T = which(mmref_temp=="T")
      loc_A = which(mmref_temp=="A")
      snp_C2T = which(mmref_temp=="C"& mm =="T")
      snp_G2A = which(mmref_temp=="G"& mm =="A")
      
      APO_TC2TT = intersect(loc_T+1,snp_C2T)
      APO_GA2AA = intersect(loc_A-1,snp_G2A)
      mm_SNP = which(mm != mmref_temp & mm != "-" & mmref_temp != "-")
      mm_SNP = setdiff(mm_SNP,APO_TC2TT)
      mm_SNP = setdiff(mm_SNP,APO_GA2AA)
      name = gsub(">","",name)
      
      TC2TT_count = length(APO_TC2TT)
      GA2AA_count = length(APO_GA2AA)
      other_snp_count = length(mm_SNP)
      all_snp_count = length(c(mm_SNP,APO_TC2TT,APO_GA2AA))
      
      if (length(mm_SNP) > 0) {
        write.table(cbind("other_SNP", mm_SNP, mmref[mm_SNP], mm[mm_SNP],other_snp_count,all_snp_count,other_snp_count/all_snp_count,other_snp_count/all_snp_count,name), paste0(filename,"_vs_",key_pattern,"_APO_info.txt"), col.names = F, 
                    row.names = F, sep = "\t", quote = F, append = T)
      }
      if (length(APO_TC2TT) > 0) {
        write.table(cbind("TC2TT", APO_TC2TT, "TC", "TT",TC2TT_count,all_snp_count,TC2TT_count/all_snp_count,(TC2TT_count+GA2AA_count)/all_snp_count, name), paste0(filename,"_vs_",key_pattern,"_APO_info.txt"), col.names = F, 
                    row.names = F, sep = "\t", quote = F, append = T)
      }
      if (length(APO_GA2AA) > 0) {
        write.table(cbind("GA2AA", APO_GA2AA, "GA", "AA",GA2AA_count,all_snp_count,GA2AA_count/all_snp_count,(TC2TT_count+GA2AA_count)/all_snp_count, name), paste0(filename,"_vs_",key_pattern,"_APO_info.txt"), col.names = F, 
                    row.names = F, sep = "\t", quote = F, append = T)
      }
      #}
    }
  }
  close(con)
}
