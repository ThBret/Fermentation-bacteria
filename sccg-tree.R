library(ape)
library(phangorn)
library(graphics)
library(TreeTools)

# Setup
Ferments = c("Acetobacter_aceti","Acetobacter_ascendens","Acetobacter_cerevisiae","Acetobacter_conturbans","Acetobacter_estunensis","Acetobacter_fabarum","Acetobacter_fallax","Acetobacter_farinalis","Acetobacter_garciniae","Acetobacter_ghanensis","Acetobacter_lambici","Acetobacter_malorum","Acetobacter_musti","Acetobacter_nitrogenifigens",
             "Acetobacter_oeni","Acetobacter_orleanensis","Acetobacter_oryzoeni","Acetobacter_pasteurianus","Acetobacter_senegalensis","Acetobacter_sicerae","Acetobacter_syzygii","Acidomonas_methanolica","Gluconacetobacter_entanii","Gluconobacter_japonicus","Gluconobacter_potus","Gluconobacter_vitians","Komagataeibacter_europaeus",
             "Komagataeibacter_intermedius","Komagataeibacter_kakiaceti","Komagataeibacter_melomenusus","Komagataeibacter_nataicola","Komagataeibacter_oboediens","Komagataeibacter_rhaeticus","Novacetimonas_hansenii","Novacetimonas_maltaceti","Acetobacter_sp","Komagataeibacter_sp","Gluconobacter_cerevisiae")

Plants = c("Acetobacter_cibinongensis","Acetobacter_indonesiensis","Acetobacter_lovaniensis","Acetobacter_okinawensis","Acetobacter_orientalis","Acetobacter_papayae","Acetobacter_peroxydans","Acetobacter_sacchari","Acetobacter_suratthaniensis","Acetobacter_thailandicus","Acetobacter_vaccinii","Ameyamaea_chiangmaiensis",
           "Asaia_astilbis","Asaia_bogorensis","Asaia_krungthepensis","Asaia_lannensis","Asaia_platycodi","Asaia_siamensis","Asaia_spathodeae","Endobacter_medicaginis","Gluconacetobacter_aggeris","Gluconacetobacter_asukensis","Gluconacetobacter_azotocaptans","Gluconacetobacter_diazotrophicus","Gluconacetobacter_dulcium",
           "Gluconacetobacter_johannae","Gluconacetobacter_liquefaciens","Gluconacetobacter_sacchari","Gluconacetobacter_takamatsuzukensis","Gluconacetobacter_tumulicola","Gluconacetobacter_tumulisoli","Gluconobacter_aidae","Gluconobacter_albidus","Gluconobacter_cadivus","Gluconobacter_frateurii","Gluconobacter_kanchanaburiensis",
           "Gluconobacter_kondonii","Gluconobacter_roseus","Gluconobacter_sphaericus","Gluconobacter_thailandicus","Gluconobacter_wancherniae","Komagataeibacter_diospyri","Komagataeibacter_medellinensis","Komagataeibacter_saccharivorans","Komagataeibacter_sucrofermentans","Komagataeibacter_swingsii","Kozakia_baliensis","Lichenicoccus_roseus",
           "Lichenicola_cladoniae","Neoasaia_chiangmaiensis","Neokomagataea_anthophila","Neokomagataea_tanensis","Neokomagataea_thailandica","Nguyenibacter_vanlangensis","Novacetimonas_cocois","Saccharibacter_floricola","Swaminathania_salitolerans","Swingsia_samuiensis","Tanticharoenia_sakaeratensis","Lichenicoccus_sp","Nguyenibacter_sp")

Insects = c("Acetobacteraceae_bacterium_ESL0697","Aristophania_vespae","Bombella_apis","Bombella_dulcis","Bombella_favorum","Bombella_intestini","Bombella_mellum","Bombella_pluederhausensis","Bombella_pollinis","Bombella_saccharophila","Bombella_sp","Commensalibacter_communis","Commensalibacter_melissae","Commensalibacter_papalotli",
            "Entomobacter_blattae","Formicincola_oecophyllae","Oecophyllibacter_saccharovorans","Parasaccharibacter_apium","Saccharibacter_sp","Commensalibacter_sp")

Fly = c("Acetobacter_oryzifermentans","Acetobacter_persici","Acetobacter_pomorum","Acetobacter_tropicalis","Commensalibacter_intestini","Gluconobacter_cerinus","Gluconobacter_morbifer","Gluconobacter_sp","Asaia_sp")

Acidophilic = c("Roseomonas_mucosa","Granulibacter_bethesdensis")

Ferments_Industrials = c("Gluconobacter_oxydans","Komagataeibacter_xylinus","Novacetimonas_pomaceti")

Sourdough = c("SRR_Acetobacter_orientalis","SRR_Acetobacter_cerevisiae","SRR_Acetobacter_fabarum_X2","SRR_Acetobacter_ghanensis_X2","SRR_Acetobacter_malorum_A1","SRR_Acetobacter_malorum_A3","SRR_Acetobacter_malorum_X1",
                    "SRR_Acetobacter_malorum_X2","SRR_Acetobacter_oryzifermentans_X1","SRR_Acetobacter_oryzifermentans_X2","SRR_Acetobacter_oryzifermentans_X3","SRR_Acetobacter_oryzifermentans_X4","SRR_Acetobacter_oryzifermentans_X5","SRR_Acetobacter_oryzoeni","SRR_Acetobacter_pasteurianus_X1",
                    "SRR_Acetobacter_pasteurianus_X2","SRR_Acetobacter_senegalensis_X1","SRR_Acetobacter_senegalensis_X2","SRR_Acetobacter_sp_X1","SRR_Acetobacter_sp_X2","SRR_Acetobacter_sp_X3","SRR_Acetobacter_syzygii","SRR_Acetobacter_tropicalis","SRR_Gluconobacter_oxydans","SRR_Gluconobacter_oxydans_B",
                    "SRR_Acetobacter_malorum_A2","SRR_Acetobacter_fabarum_X1","SRR_Acetobacter_ghanensis_X1","SRR_Acetobacter_okinawensis")

# PLOT SCCG TREE
pdf("/Users/thibaultbret/sccg-tree.pdf", width = 32, height = 28)  # Adjust width and height as needed (40 x 45 with SRR; 32 x 38 without)

plot_indiv_tree <- function(input_tree, id, bootstrap = FALSE, SRR = FALSE, time_scale = FALSE, plot = TRUE, aligned = FALSE){
  tree_file <- read.tree(input_tree)
  if(aligned==TRUE){tree_file <- drop.tip(tree_file, c("Acetobacteraceae_bacterium_ESL0697","Aristophania_vespae","Bombella_saccharophila","Bombella_pollinis","Bombella_sp","Bombella_apis","Parasaccharibacter_apium","Bombella_dulcis"))}
  label_colours <- rep("black", length(tree_file$tip.label))
  label_colours[tree_file$tip.label %in% Ferments] <- "#9933CC"
  label_colours[tree_file$tip.label %in% Insects] <- "#E6A401"
  label_colours[tree_file$tip.label %in% Plants] <- "#73C221"
  label_colours[tree_file$tip.label %in% Fly] <- "#CCCCCC"
  label_colours[tree_file$tip.label %in% Ferments_Industrials] <- "#FF00FF"
  label_colours[tree_file$tip.label %in% Sourdough] <- "#0356FC"
  label_colours[tree_file$tip.label %in% Acidophilic] <- "#299596"
  
  # Time scale
  if(time_scale==TRUE){
    #getMRCA to figure out node numbers
    getMRCA(tree_file, c("Aristophania_vespae","Parasaccharibacter_apium")) # node 167

    # Calibrate the time scale based on no specified data
    #calib <- makeChronosCalib(tree_file, node = "root", age.max = 1)
    # Calibrate based on a split of node 167 estimated to have taken place 80 million years ago
    calib <- makeChronosCalib(tree_file, node = 167, age.max = 80)
    time_tree <- chronos(tree_file, lambda = 1, calibration = calib, model = "discrete")
    # Plot the tree without bootstrap values
    plot(time_tree, edge.color = "black", tip.color = label_colours, main = "Time-calibrated single-Copy Core Genes Tree", cex.main = 3.5, cex = 1.5, show.tip.label = TRUE)
    # Add a time scale
    axisPhylo(cex = 2)
    graphics::mtext("Time (myrs)", side = 1, line = 3, cex = 2, at = max(get("last_plot.phylo",envir = .PlotPhyloEnv)$xx) * 0.5)
  } else {
    # Plot the tree with bootstrap values
    if(plot==TRUE){
      if(bootstrap==TRUE){
        plotBS(tree_file, edge.color = "black", tip.color = label_colours, main = "Single-Copy Core Genes Tree", cex.main = 3.5, cex = 1.5)
        } else {
          if(aligned==FALSE){
            plot(tree_file, edge.color = "black", tip.color = label_colours, main = "Single-Copy Core Genes Tree", cex.main = 3.5, cex = 1.5)
            node_colours <- ifelse(tree_file$node.label < 33.33, "red", NA)
            nodelabels(text = character(length(node_colours)), bg = node_colours, frame = "none", pch = ifelse(!is.na(node_colours), 21, NA), cex = 2)
          } else {
            plot(tree_file, direction = "leftwards", align.tip.label = TRUE, x.lim = 15, edge.color = "black", tip.color = label_colours, ps = 12)
            }
          }
    } else {
      return(label_colours)
    }
  }
  
  # Legend
  if(aligned==TRUE){
    #No legend
  } else {
  if(SRR==TRUE){
    legend("topright",
           legend = c("Ferment", "Plant", "Insect","Fly","Ferment/Industrial","Sourdough"), #Name of groups
           fill = c("#9933CC","#73C221","#E6A401","#CCCCCC","#FF00FF","#0356FC"), # Colour of the squares
           border = "black", # Colour of the border of the squares
           cex = 3, #sets legend size
           xpd = FALSE) #places outside plot area
  }else{
    legend("topright",
           legend = c("Ferment", "Plant", "Insect","Fly","Ferment/Industrial"), #Name of groups
           fill = c("#9933CC","#73C221","#E6A401","#CCCCCC","#FF00FF"), # Colour of the squares
           border = "black", # Colour of the border of the squares
           cex = 3, #sets legend size
           xpd = FALSE) #places outside plot area
    }
  }
}

plot_indiv_tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile", SRR = FALSE, time_scale = FALSE,  aligned = TRUE)
#plot_indiv_tree("/Users/thibaultbret/Documents/Work/concord_withSRR.treefile", bootstrap = TRUE, SRR = FALSE, time_scale = FALSE)

#nodelabels(cex = 0.4)

dev.off()



########################################################################################################################
## PAIRWISE COMPARISONS

# Komagateibacter vs Gluconacetobacter
library(TreeTools)
tree <- read.tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile")
tree <- Preorder(tree)
getMRCA(tree, c("Nguyenibacter_sp","Novacetimonas_pomaceti")) # node 139

write.tree(Subtree(tree, 139), "/Users/thibaultbret/Komagateibacter-Gluconacetobacter-tree.treefile")

# Acetobacter clades vs Gluconobacter clades + Neokomagatea clade
getMRCA(tree, c("Acetobacter_tropicalis","Gluconobacter_wancherniae")) # node 170
subtree <- Subtree(tree, 140)
subtree <- drop.tip(subtree, c("Bombella_apis","Bombella_dulcis","Bombella_favorum","Bombella_intestini","Bombella_mellum","Bombella_pluederhausensis","Bombella_pollinis",
                    "Bombella_saccharophila","Bombella_sp","Formicincola_oecophyllae","Oecophyllibacter_saccharovorans","Aristophania_vespae","Acetobacteraceae_bacterium_ESL0697",
                    "Parasaccharibacter_apium","Saccharibacter_floricola","Saccharibacter_sp","Asaia_sp","Kozakia_baliensis","Swaminathania_salitolerans","Asaia_astilbis",
                    "Asaia_bogorensis","Asaia_krungthepensis","Asaia_lannensis","Asaia_platycodi","Asaia_siamensis","Asaia_spathodeae","Acidomonas_methanolica","Neoasaia_chiangmaiensis",
                    "Ameyamaea_chiangmaiensis","Tanticharoenia_sakaeratensis","Parasaccharibacter_apium"))

write.tree(subtree, "/Users/thibaultbret/Aceto_vs_GluconoNeoko-tree.treefile")

#
getMRCA(tree, c("Entomobacter_blattae","Bombella_mellum")) # node 138
subtree <- Subtree(tree, 138)
subtree <- drop.tip(subtree, c("Gluconobacter_kanchanaburiensis","Gluconobacter_cerevisiae","Gluconobacter_kondonii","Gluconobacter_cadivus","Gluconobacter_albidus","Gluconobacter_sphaericus","Gluconobacter_aidae","Gluconobacter_sp","Gluconobacter_oxydans","Gluconobacter_potus","Gluconobacter_roseus","Gluconobacter_vitians",
                               "Gluconobacter_morbifer","Gluconobacter_thailandicus","Gluconobacter_frateurii","Gluconobacter_japonicus","Gluconobacter_cerinus","Neokomagataea_tanensis","Neokomagataea_thailandica","Neokomagataea_anthophila","Swingsia_samuiensis","Gluconobacter_wancherniae","Asaia_bogorensis","Asaia_krungthepensis",
                               "Asaia_sp","Asaia_lannensis","Asaia_spathodeae","Asaia_siamensis","Asaia_platycodi","Asaia_astilbis","Swaminathania_salitolerans","Kozakia_baliensis","Neoasaia_chiangmaiensis","Acidomonas_methanolica","Tanticharoenia_sakaeratensis","Ameyamaea_chiangmaiensis","Acetobacter_oryzoeni","Acetobacter_pomorum",
                               "Acetobacter_oryzifermentans","Acetobacter_ascendens","Acetobacter_pasteurianus","Acetobacter_vaccinii","Acetobacter_papayae","Acetobacter_suratthaniensis","Acetobacter_peroxydans","Acetobacter_garciniae","Acetobacter_okinawensis","Acetobacter_lambici","Acetobacter_lovaniensis","Acetobacter_fabarum","Acetobacter_syzygii",
                               "Acetobacter_ghanensis","Acetobacter_cibinongensis","Acetobacter_orientalis","Acetobacter_thailandicus","Acetobacter_indonesiensis","Acetobacter_tropicalis","Acetobacter_senegalensis","Acetobacter_cerevisiae","Acetobacter_malorum","Acetobacter_orleanensis","Acetobacter_persici","Acetobacter_farinalis",
                               "Acetobacter_musti","Acetobacter_fallax","Acetobacter_oeni","Acetobacter_conturbans","Acetobacter_aceti","Acetobacter_sicerae","Acetobacter_sp","Acetobacter_estunensis","Acetobacter_nitrogenifigens","Acetobacter_sacchari","Gluconacetobacter_azotocaptans","Gluconacetobacter_johannae","Gluconacetobacter_tumulisoli","Gluconacetobacter_diazotrophicus",
                               "Nguyenibacter_sp","Nguyenibacter_vanlangensis","Gluconacetobacter_asukensis","Gluconacetobacter_tumulicola","Gluconacetobacter_aggeris","Gluconacetobacter_takamatsuzukensis","Gluconacetobacter_liquefaciens","Gluconacetobacter_dulcium","Gluconacetobacter_sacchari","Komagataeibacter_rhaeticus","Komagataeibacter_medellinensis","Komagataeibacter_swingsii",
                               "Komagataeibacter_europaeus","Komagataeibacter_diospyri","Komagataeibacter_intermedius","Komagataeibacter_oboediens","Komagataeibacter_melomenusus","Komagataeibacter_xylinus","Komagataeibacter_sucrofermentans","Komagataeibacter_nataicola","Komagataeibacter_saccharivorans","Komagataeibacter_sp","Komagataeibacter_kakiaceti","Gluconacetobacter_entanii",
                               "Novacetimonas_maltaceti","Novacetimonas_hansenii","Novacetimonas_pomaceti","Novacetimonas_cocois"))

write.tree(subtree, "/Users/thibaultbret/Bombe_vs_Commens-tree.treefile")



########################################################################################################################
# ANCESTRAL METABOLIC RECONSTRUCTION TREES

library(phytools)
tree <- read.tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile")
mean_presence <- read.table("/Users/thibaultbret/Mean-expression-rates.txt", header = TRUE, row.names = 1)
mean_presence$normalised_presence = (mean_presence$Total-min(mean_presence$Total))/(max(mean_presence$Total)-min(mean_presence$Total))
amino_acid_expression_rates <- read.table("/Users/thibaultbret/Amino-acid-expression-rates.txt", header = TRUE, row.names = 1)
amino_acid_expression_rates$normalised_presence = (amino_acid_expression_rates$Amino_acid_metabolism-min(amino_acid_expression_rates$Amino_acid_metabolism))/(max(amino_acid_expression_rates$Amino_acid_metabolism)-min(amino_acid_expression_rates$Amino_acid_metabolism))
carbohydrate_expression_rates <- read.table("/Users/thibaultbret/Carbohydrate-expression-rates.txt", header = TRUE, row.names = 1)
carbohydrate_expression_rates$normalised_presence = (carbohydrate_expression_rates$Carbohydrate_metabolism-min(carbohydrate_expression_rates$Carbohydrate_metabolism))/(max(carbohydrate_expression_rates$Carbohydrate_metabolism)-min(carbohydrate_expression_rates$Carbohydrate_metabolism))

svl <- as.matrix(mean_presence[3])[,1]
#svl <- as.matrix(amino_acid_expression_rates[3])[,1]
#svl <- as.matrix(carbohydrate_expression_rates[3])[,1]

fit <- fastAnc(tree, svl, vars=TRUE, CI=TRUE)

## projection of the reconstruction onto the edges of the tree
pdf("/Users/thibaultbret/ancestral-metabolic-reconstruction24.pdf", width = 18, height = 24)  # Adjust width and height as needed
obj <- contMap(tree, svl, plot = FALSE, res = 150)
col_scale <- c("#f0f921","#fca636","#e16462","#b12a90","#6a00a8","#0d0887")
obj <- setMap(obj, colors = col_scale)
plot(obj, type="fan", legend = max(nodeHeights(tree)), fsize=c(0.7,0.9), offset = 3, leg.txt="Normalised overall presence estimate of KEGG amino acid metabolism pathways")
dev.off()



########################################################################################################################
## TIME TREE WITH GEOLOGICAL SCALE
library(strap)
library(MCMCtreeR)
calibrated_tree <- readMCMCtree("/Users/thibaultbret/Documents/Work/calibrated-mcmc-tree.tre")
geoscalePhylo(calibrated_tree)

#
tree <- read.tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile")
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Annotation rows
label_colours <- c("#299596","#299596","#73C221","#9933CC","#73C221","#73C221","#73C221","#73C221","#73C221","#CCCCCC","#FF00FF","#9933CC","#73C221","#9933CC","#CCCCCC","#73C221","#73C221","#9933CC","#CCCCCC","#73C221","#73C221","#73C221","#73C221","#73C221","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#73C221","#E6A401","#E6A401","#E6A401","#E6A401","#E6A401","#73C221","#73C221","#CCCCCC","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#9933CC","#73C221","#73C221","#9933CC","#CCCCCC","#CCCCCC","#9933CC","#9933CC","#73C221","#73C221","#73C221","#73C221","#9933CC","#73C221","#9933CC","#73C221","#9933CC","#9933CC","#9933CC","#73C221","#73C221","#73C221","#73C221","#CCCCCC","#9933CC","#9933CC","#9933CC","#9933CC","#CCCCCC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#73C221","#9933CC","#73C221","#73C221","#9933CC","#73C221","#9933CC","#9933CC","#9933CC","#FF00FF","#73C221","#9933CC","#73C221","#9933CC","#9933CC","#9933CC","#9933CC","#9933CC","#FF00FF","#73C221","#E6A401","#E6A401","#CCCCCC","#E6A401","#E6A401","#E6A401","#73C221","#73C221","#73C221","#73C221","#000000")

pdf("/Users/thibaultbret/geological-tree.pdf", width = 15, height = 15)  # Adjust width and height as needed
MCMC.tree.plot(calibrated_tree, cex.tips = .6, time.correction = 100, scale.res = c("Eon","Period", "Epoch"), 
            plot.type = "phylogram", cex.age = 0.8, cex.labels = 0.7, relative.height = 0.08, col.tree = "grey40", label.offset = 1, 
            node.method = "bar", grey.bars = FALSE, col.age = "#10179c50", lwd.bar = 8, no.margin = TRUE, tip.color = label_colours)
dev.off()

#just node 168
pdf("/Users/thibaultbret/geological-tree-highlighted-split.pdf", width = 15, height = 15)  # Adjust width and height as needed
MCMC.tree.plot(calibrated_tree, cex.tips = .6, time.correction = 100, scale.res = c("Eon","Period", "Epoch"), 
               plot.type = "phylogram", cex.age = 0.8, cex.labels = 0.7, relative.height = 0.08, col.tree = "grey40", label.offset = 1, 
               node.method = "bar", col.age = "#ff000040", lwd.bar = 30, all.nodes = 168, no.margin = TRUE)
dev.off()


#cladogram
pdf("/Users/thibaultbret/geological-cladogram.pdf", width = 15, height = 15)  # Adjust width and height as needed
MCMC.tree.plot(calibrated_tree, analysis.type = "MCMCtree", cex.tips = 0.6, label.offset = 1,
               time.correction = 100, plot.type = "cladogram", lwd.bar = 2, 
               scale.res = c("Eon", "Period","Epoch"), node.method = "node.length", 
               col.age = "#008b0080", no.margin = TRUE, cex.labels = .7)
dev.off()



########################################################################################################################
## PHYLOGENETIC SIGNAL
library(phytools)
tree <- read.tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile")
extradata <- read.delim("/Users/thibaultbret/extradata.tsv", row.names = 1) 
## extract discrete characters
log_genome_sizes <- setNames(log(extradata$Size),
                        rownames(extradata))
log_GC_contents <-  setNames(log(extradata$GC.content),
                          rownames(extradata))
#cont
genome_sizes_K <- phylosig(tree, log_genome_sizes, test = TRUE)
genome_sizes_lambda <- phylosig(tree, log_genome_sizes, method = "lambda", test = TRUE)
GC_contents_K <- phylosig(tree, log_GC_contents, test = TRUE)
GC_contents_lambda <- phylosig(tree, log_GC_contents, method = "lambda", test = TRUE)


# plot results
pdf("/Users/thibaultbret/phylogenetic-signals.pdf", width = 10, height = 10)  # Adjust width and height as needed
par(mfrow=c(2,2),cex=0.9)
plot(genome_sizes_K,las=1,cex.axis=0.9)
mtext("a)",adj=0,line=1)
plot(genome_sizes_lambda,bty="n",las=1,cex.axis=0.9,
     xlim=c(0,1.1))
mtext("b)",adj=0,line=1)
plot(GC_contents_K,las=1,cex.axis=0.9)
mtext("c)",adj=0,line=1)
plot(GC_contents_lambda,bty="n",las=1,cex.axis=0.9,
     xlim=c(0,1.1))
mtext("d)",adj=0,line=1)
dev.off()

# Genome size tree
genome_sizes <- setNames(extradata$Size / 1000000,
                         rownames(extradata))
genomesize.widthMap <- edge.widthMap(tree, genome_sizes)
pdf("/Users/thibaultbret/genome-size-tree.pdf", width = 18, height = 24)  # Adjust width and height as needed
plot(genomesize.widthMap, color = palette()[2], no.margin=TRUE, fsize = .8)
dev.off()

# GC content tree
gc_contents <- setNames(extradata$GC.content,
                         rownames(extradata))
gc_contents.widthMap <- edge.widthMap(tree, gc_contents)
pdf("/Users/thibaultbret/gc-content-tree.pdf", width = 18, height = 24)  # Adjust width and height as needed
plot(gc_contents.widthMap, color = palette()[4], no.margin=TRUE, fsize = .8)
dev.off()

# Hidden transition states
# remove acidophilic & human disease to restrict genomes to 4 isolation sources
extradata_4cat <- extradata[!(extradata$Isolation.source %in% c('Acidophilic','Human disease')),]
tree_4cat <- drop.tip(tree, c("Acidocella_aromatica","Roseomonas_mucosa","Granulibacter_bethesdensis"))
isolation_sources <- setNames(extradata_4cat$Isolation.source,
                              rownames(extradata_4cat))

# MK model
fitMK <- fitMk(tree_4cat, isolation_sources, model = "ARD")

# plot results
pdf("/Users/thibaultbret/that.pdf", width = 24, height = 18)  # Adjust width and height as needed

fitMK_asr <- ancr(fitMK, tips = TRUE)
cols<-setNames(c("#9933CC","#CCCCCC","#E6A401",
                 "#73C221"),colnames(fitMK_asr$ace))
plot(fitMK_asr,legend=FALSE,
     args.plotTree = list(type="arc",offset = 5, arc_height= .5),
     args.nodelabels = list(piecol = cols, cex = 0.2),
     args.tiplabels = list(cex = 0.05))

#
xy<-plot(as.Qmatrix(fitMK),show.zeros = FALSE,
    ylim=c(-1.2,3),xlim=c(-2,2),spacer = 0.2,
         text = TRUE, width = TRUE, color = TRUE, max.lwd = 5)
nulo<-mapply(plotrix::draw.circle,xy$x,xy$y,col=cols,
             MoreArgs=list(radius=0.15))
#text(0,0.1,"fitted MK model",adj=0.5, cex = 2)
text(xy$x,xy$y,xy$states,cex= 1.5,col="white")

dev.off()
