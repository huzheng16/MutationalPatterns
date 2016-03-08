# -------------- MAKE METADATA FILE -----------------------

sample_info = read.table("~/Dropbox/francis/mutational_patterns_ASC/data/sample_info.txt", sep = "\t", header =T)
donor_info =  read.table("~/Dropbox/francis/mutational_patterns_ASC/data/donor_info.txt", sep = "\t", header =T)

# INTESTINE
surveyed_files_intestine = list.files("~/Documents/Organoids/ASC_multi_tissues/data/surveyed/intestine/", full.names = T)
vcf_files_intestine = list.files("~/Documents/Organoids/ASC_multi_tissues/data/vcf_filtered/intestine/", full.names = T)
# LIVER
surveyed_files_liver = list.files("~/Documents/Organoids/ASC_multi_tissues/data/surveyed/liver/", full.names = T)
vcf_files_liver = list.files("~/Documents/Organoids/ASC_multi_tissues/data/vcf_filtered/liver/", full.names = T)
# COLON
surveyed_files_colon = list.files("~/Documents/Organoids/ASC_multi_tissues/data/surveyed/colon/", full.names = T)
vcf_files_colon = list.files("~/Documents/Organoids/ASC_multi_tissues/data/vcf_filtered/colon/", full.names = T)

# vcf files
vcf_files = c(vcf_files_colon, vcf_files_intestine, vcf_files_liver)
samples = unlist(lapply(vcf_files, function(x) extract_sample_name(x)))
vcf_file_info = data.frame(vcf_file=vcf_files, Sample=samples)
# surveyed bed files
surveyed_files = c(surveyed_files_colon, surveyed_files_intestine, surveyed_files_liver)
samples_surveyed = unlist(lapply(surveyed_files, function(x) extract_sample_name(x)))
surveyed_file_info = data.frame(surveyed_file=surveyed_files, Sample=samples_surveyed)
file_info = merge(vcf_file_info, surveyed_file_info)
# add files to sample info
sample_info = merge(sample_info, file_info)
info = merge(donor_info, sample_info, by = c("Donor", "Donor"))

write.table(info, file ="~/Dropbox/francis/mutational_patterns_ASC/data/info.txt", sep = "\t", quote = F, col.names = T, row.names=F)

# --------------------- READ DATA IN MutPat OBJECT ----------------------------------
vcf_files = as.character(info$vcf_file)
vcf = read_vcf_list(vcf_files)
# surveyed = bed_to_granges_list(as.character(info$surveyed_files)) 
MutPat_object = list(vcf = vcf, type = info$Tissue, individual = info$Donor, age = info$Age)


type_colors = c("#E6204E", "#F7B75B", "#B78ABE")
sub_colors = c("#DBD7C8", "#B2D39C", "#71C1BA", "#2DAFCE", "#2476B2", "#737E93")

