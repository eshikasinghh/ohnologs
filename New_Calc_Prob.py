import sys

ref_genome_name = sys.argv[1]
ref_genome_file = sys.argv[2]
sub_genome_name = sys.argv[3]
sub_genome_file = sys.argv[4]
orthology_file = sys.argv[5]
reference_window = int(sys.argv[6]) #converts to integer b/c should be a number
polyploid_window = int(sys.argv[7])  
min_orthologs = int(sys.argv[8])   
shrink = sys.argv[9]
out_files_dir = sys.argv[10]

ortholog_relation_file = f"{ref_genome_name}-{sub_genome_name}_Ortholog_Relations.txt"
filtered_blocks_file = f"{ref_genome_name}-{sub_genome_name}_FilteredBlocks_Window-P{polyploid_window}R{reference_window}_Orthologs-{min_orthologs}.txt"
pp_gene_count = f"{sub_genome_name}_Chromosomes-Genes.txt"
og_gene_count = f"{ref_genome_name}_Chromosomes-Genes.txt"
blocks_with_probability = f"{ref_genome_name}-{sub_genome_name}_FilteredBlocksWithProbability_Window-P{polyploid_window}R{reference_window}_Orthologs-{min_orthologs}.txt"
anchor_file_with_p = f"{ref_genome_name}-{sub_genome_name}_AnchorsWithProbability_Window-P{polyploid_window}R{reference_window}_Orthologs-{min_orthologs}.txt"
#-------------------------------------------------------------------------------------------------------------------
# In case if the window is symmatric do the following
def getSymmWindowsWithoutAnyOrtholog(totalGenes, wdSize, genes):
    # Sort the genes list
    genes = sorted(genes)
    
    #print(f"total genes = {totalGenes}\twindow = {wdSize}\tTotal orth on this chr = ", end="")
    # print(len(genes), "\tOrthologs = ")
    
    # Calculate the segments between genes from the start to the end of the chromosome
    # Append start position = 0 and end position = totalGenes + 1
    # to the list if the first and last elements are not already 1 or the length of the gene.
    genes.append(totalGenes + 1)
    
    # print('|'.join(map(str, genes)))
    
    start = 0  # first position of gene
    totalPossibleWindows = 0  # total possible windows without any ortholog of this gene
    
    for i in range(len(genes)):
        
        segment = genes[i] - start - 1
        # print(f"Start = {start}\tGenePos = {genes[i]}\tSegment = {segment}")
        
        if start == 0 or genes[i] == (totalGenes + 1):  # First or last segment
            # Window will shrink and be half for symmetric window
            if segment > (wdSize / 2):
                totalPossibleWindows += (segment - (wdSize / 2))
                # print(f"* Windows = {totalPossibleWindows}")
        
        elif start != 0 and genes[i] != (totalGenes + 1):  # Middle segment
            # The window will not shrink
            if segment > wdSize:
                totalPossibleWindows += (segment - wdSize)
                # print(f"-> Windows = {totalPossibleWindows}")

        start = genes[i]

    # print(f"Total Windows for this gene = {totalPossibleWindows}")

    probabilityOf_NONE_OfGenesWithinWindow = totalPossibleWindows / totalGenes
    probabilityOf_ONE_OfGeneWithinWindow = 1 - probabilityOf_NONE_OfGenesWithinWindow

    # print(f"{probabilityOf_NONE_OfGenesWithinWindow}\t{probabilityOf_ONE_OfGeneWithinWindow}")

    return probabilityOf_ONE_OfGeneWithinWindow

# In case of Asymmatric window, do the following
def getAsymmWindowsWithoutAnyOrtholog(totalGenes, wdSize, genes):
    # Sort the genes list
    genes = sorted(genes)
    totalGenes = int(totalGenes)
    
    # print(f"total genes = {totalGenes}\twindow = {wdSize}\tTotal orth on this chr = ", end="")
    # print(len(genes), "\tOrthologs = ")
    
    # Calculate the segments between genes from the start of the chromosome to the end
    # Append the end position = totalGenes + 1 to the list
    genes.append(totalGenes + 1)
    
    # print(">")
    # print('|'.join(map(str, genes)))
    
    start = 0  # first position of gene
    totalPossibleWindows = 0  # total possible windows WITHOUT ANY ORTHOLOG of this gene
    
    for i in range(len(genes)):
        segment = genes[i] - start - 1
        # print(f"Start = {start}\tGenePos = {genes[i]}\tSegment = {segment}")

        if segment > wdSize:  # For asymmetric windows, the window size is always the same, so no shrinking
            totalPossibleWindows += (segment - wdSize)
            # print(f"Windows = {totalPossibleWindows}")
        
        start = genes[i]
    
    # print(f"Total Windows for this gene = {totalPossibleWindows}")
    
    # print(f"{probabilityOf_NONE_OfGenesWithinWindow}\t{probabilityOf_ONE_OfGeneWithinWindow}")
    
    return totalPossibleWindows, (totalGenes - wdSize)

# The same subrouteins to get the position of the asymmatric and symmatric window to get total og genes in the window having orthologs anywhere on the genome
def getAsymmetricWindowPosition(pos, N, wd):
    # pos: position
    # N: total number of genes
    # wd: this is already a half window
    
    leftWd = pos - wd
    rightWd = pos + wd
    
    # If the left window is out of bounds, adjust it and the right window
    if leftWd <= 0:
        lOffset = abs(leftWd) + 1
        rightWd += lOffset
        leftWd = 1
    
    # Fix the right window if it exceeds the chromosome boundary
    if rightWd > N:
        rOffset = rightWd - N
        leftWd -= rOffset
        rightWd = N
    
    # Ensure the windows are within the chromosome boundaries
    leftWd = max(1, leftWd)  # Ensure left boundary is at least 1
    rightWd = min(N, rightWd)  # Ensure right boundary does not exceed N
    
    return int(leftWd), int(rightWd)  # Ensure the return values are integers

def getSymmetricWindowPosition(pos, N, wd):
    # pos: position
    # N: total number of genes
    # wd: half window size
    
    # Get the window boundaries simply by adding the half window to both sides of the gene
    # Example: if position is at 755 and window size is 200, the window will be 655 to 855 (201 genes including 655th and 855th genes)
    leftWd = pos - wd
    rightWd = pos + wd
    
    # If the window boundaries exceed the segment length, reset them to segment boundaries
    if leftWd <= 0:
        leftWd = 1
    if rightWd > N:
        rightWd = N
    
    return int(leftWd), int(rightWd)
#_______________
print("Calculating probability of random occurrence...")

# open outfiles
try:
    prob_file = open(f"{out_files_dir}/{blocks_with_probability}", "w")
except IOError as e:
    print(f"Could not open file {out_files_dir}/{blocks_with_probability} for writing: {e}")
    sys.exit(1)

try:
    ancprob_file = open(f"{out_files_dir}/{anchor_file_with_p}", "w")
except IOError as e:
    print(f"Could not open file {out_files_dir}/{anchor_file_with_p} for writing: {e}")
    sys.exit(1)

#----- get the sum for all Ks that herve wants --------   NOT USING RIGHT NOW
#open K1, '>Sum of P1 for all k.txt' or die $!;
#open K2, '>Sum of P2 for all k.txt' or die $!;


# read ortholog relation file
try:
    with open(f"{out_files_dir}/{ortholog_relation_file}", "r") as or_file:
        orf = or_file.readlines()
except IOError as e:
    print(f"Could not open file {out_files_dir}/{ortholog_relation_file} for reading: {e}")
    sys.exit(1)
    
orf = orf[1:]

# and make the superhash. I need it to get all the orthologs not just ones in the blocks

super_dict = {}

# Process each line in ortholog relation file
for each in orf:
    line = each.strip().split("\t")
    line = [elem.replace("\n", "") for elem in line]

    if all(line[:4]):
        key1, key2, key3, key4 = line[:4]
       
        if key1 not in super_dict:
            super_dict[key1] = {}
        if key2 not in super_dict[key1]:
            super_dict[key1][key2] = {}
        if key3 not in super_dict[key1][key2]:
            super_dict[key1][key2][key3] = {}
        super_dict[key1][key2][key3][key4] = ''

# Process each line in ortholog relation file



# Read WGD chromosome and gene count file
try:
    with open (f"{out_files_dir}/{pp_gene_count}", "r") as pgcfile:  #with statement automatically closes file
        pgc = pgcfile.readlines()
except IOError as e:
    print(f"Cannot open file: {e}")

PolyploidGeneCount = {}
for each in pgc:
    line = each.strip().split("\t")
    PolyploidGeneCount[line[0]] = line[1]

# Read outgroup chromosome and gene count file
try:
    with open(f"{out_files_dir}/{og_gene_count}", "r") as ogcfile:
        ogc = ogcfile.readlines()
except IOError as e:
    print(f"Cannot open file: {e}")

OutgroupGeneCount = {}
for each in ogc:
    line = each.strip().split("\t")

    # chromosome_name = line[0].strip()
    # if not chromosome_name.startswith('chr'):
    #     chromosome_name = 'chr' + chromosome_name

    # Store the cleaned chromosome name and gene count
    #OutgroupGeneCount[chromosome_name] = line[1].strip()

# Debugging step to ensure correct loading of data
#print("OutgroupGeneCount keys:", OutgroupGeneCount.keys())
    OutgroupGeneCount[line[0]] = line[1]
#print join "\t", keys %OutgroupGeneCount;
#print("\t".join(sorted(OutgroupGeneCount.keys())))


# Read filtered blocks 
    #python doesn't have proper equivalent to local $/ below is chatgpt
with open(f"{out_files_dir}/{filtered_blocks_file}", 'r') as fh:
    content = fh.read()

blocks = content.split('>')

if blocks[0].strip() == "":
    blocks.pop(0)

for block in blocks:
    genes_in_this_block = block.split("\n")
    if '>' in block:
        genes_in_this_block.pop()
    anchor = genes_in_this_block.pop(0) if genes_in_this_block else None
    genes_in_block_with_probability = []

#--- A few variables to store log p values for the block, and calculate final probability ----------------------#
    p_all_genes_within_window = 1           # Multiplication of Ps for orthologs in Og-PP window for entire genome -- EXCLUDING ANCHOR
    p_one_genes_within_window_chr = None    # Ps for orthologs in Og-PP window for only current chromosomes
    p_all_genes_within_window_chr = 1       # Multiplication of Ps for orthologs in Og-PP window for only current chromosoms -- EXCLUDING ANCHOR
    p_anchor_genome = None                  # P for all genes anchor FOR ENTIRE GENOME
    p_anchor_chr = None                     # P for all genes anchor FOR JUST CURRENT CHROMOSOME
    sum_log_p = 0                           # Summation of log of probability values 
    sum_log_1_minus_p = 0                   # Summation of log of 1 minus probability values
    og_genes_with_orth_in_this_wd = {}      # Hash to identify unique Og genes having orthologs in this window ACROSS ENTIRE GENOME
    og_genes_with_orth_in_this_wd_chr = {}  # Hash to identify unique Og genes having orthologs in this window ONLY ON CURRENT CHROMOSOME
    k = None                                # Total UNIQUE og genes having orthologs in current PP window -- EXCLUDING ANCHOR
    N = None                                # Total UNIQUE og genes having orthologs anywhere in PP chromosomes -- EXCLUDING ANCHOR
	
#	print "$anchor\n";  
    # print(f"{anchor}")

#uses regex patterns, chatgpt below
    import re
    anchor = anchor.strip()
    #anchor = "Anchor: chr1\t12345\tchr2\t67890\t...\t...\t...\t..."
    #pattern = r"Anchor: (.+)\t(\d+)\t(.+)\t(\d+)\t.+\t.+\t.+\t.+"
    pattern = r"Anchor: (\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.+)"

    match = re.search(pattern, anchor)

    if match:  # Extract captured groups
        og_anc_chr = match.group(1)   # First capture group
        og_anc_pos = int(match.group(2))  # Second capture group, converted to integer
        pp_anc_chr = match.group(3)   # Third capture group
        pp_anc_pos = int(match.group(4)) # Fourth capture group, converted to integer
    
    #print(f"Raw anchor line: {anchor}")
    #print(f"{og_anc_chr} {og_anc_pos} {pp_anc_chr} {pp_anc_pos}")
    #print(f"{og_anc_pos}, {OutgroupGeneCount[og_anc_chr]}, {int(reference_window / 2)}")
    
#--- Get the window boundary for outgroup ---#
    
    if shrink == 'n':
        leftWdPos, rightWdPos = getAsymmetricWindowPosition(og_anc_pos, int(OutgroupGeneCount[og_anc_chr]), reference_window / 2)
        #print(f"Shrink: 'n', Left Window Position: {leftWdPos}, Right Window Position: {rightWdPos}")
    elif shrink == 'y':
        leftWdPos, rightWdPos = getSymmetricWindowPosition(og_anc_pos, int(OutgroupGeneCount[og_anc_chr]), reference_window / 2)
        #print(f"Shrink: 'y', Left Window Position: {leftWdPos}, Right Window Position: {rightWdPos}")
    
    leftWdPos, rightWdPos = int(leftWdPos), int(rightWdPos)
    #--- For each outgroup gene position in window do the following ---#
    uniqueGenesWithOrth = 0  # This is the total number of outgroup genes in the current window having orthologs in the WGD genome

# [Previous code remains the same until the window position loop]

    for position in range(leftWdPos, rightWdPos + 1):  # For each position
        if og_anc_chr in super_dict and str(position) in super_dict[og_anc_chr]:  # If it has orthologs
            uniqueGenesWithOrth += 1
            
            # Initialize dictionary with all possible polyploid chromosomes as keys
            allPPChrPos = {ppChr: [] for ppChr in PolyploidGeneCount.keys()}
            
            # Fill in the positions for chromosomes that have orthologs
            for ppChr in super_dict[og_anc_chr][str(position)].keys():
                # Convert the dictionary keys (positions) to a sorted list
                positions = sorted([int(pos) for pos in super_dict[og_anc_chr][str(position)][ppChr].keys()])
                allPPChrPos[ppChr] = positions
            
            # Format output string similar to Perl
            chr_pos_pairs = []
            for ppChr, positions in sorted(allPPChrPos.items()):
                if positions:  # Only include chromosomes that have positions
                    pos_str = f"[{', '.join(map(str, positions))}]"
                    chr_pos_pairs.append(f"'{ppChr}' => {pos_str}")
            
            #output_str = f"{og_anc_chr} - {position} => " + "{" + ", ".join(chr_pos_pairs) + "}"
            #print(output_str)
            totalWdsWithoutAnyGenes = 0
            totalWdsonALLChromosomes = 0     #line 184

            for chromosome in allPPChrPos.keys():
    # Pass the gene count on 'chromosome' and all genes on 'chromosome'
    # If there are no genes on any chromosome, the empty list will be passed, as declared above
                totWd = 0
                wdNoGene = 0
    
                if shrink == 'y':
                    # Get the probability that this gene lies within any window by chance
                    wdNoGene, totWd = getSymmWindowsWithoutAnyOrtholog(PolyploidGeneCount[chromosome], polyploid_window, allPPChrPos[chromosome])  # Pass total genes on current chromosome, window size, and list having all gene positions on this chromosome
    
                if shrink == 'n':
                    # Get the probability that this gene lies within any window by chance
                    wdNoGene, totWd = getAsymmWindowsWithoutAnyOrtholog(PolyploidGeneCount[chromosome], polyploid_window, allPPChrPos[chromosome])  # Pass total genes on current chromosome, window size, and list having all gene positions on this chromosome
                # printing outgroup chromosome|position and position of orthologs on all polyploid chromosomes if there is any ortholog, and possible windows
				    
                # print(f"{og_anc_chr}|{position}\t{ppChr}\t|{', '.join(map(str, allPPChrPos[ppChr]))}|\t{wdNoGene}\t{totWd}")
                
                if totWd >= 0:
                    totalWdsWithoutAnyGenes += wdNoGene
                    totalWdsonALLChromosomes += totWd
                 # Exiting after print for debugging


                 # Exiting after print for debugging
# If the polyploid anchor is the current chromosome, take it in a new variable 
# to get the probability that we will find at least 1 orth ONLY ON CURRENT POLYPLOID CHROMOSOME randomly in the window
                if pp_anc_chr == chromosome:
                    #print(f"{chromosome}\t|{', '.join(map(str, allPPChrPos[chromosome]))}|\t{wdNoGene}\t{totWd}")
                   # Exiting as in the Perl script for debugging purposes
    # print(f"{chromosome}\t{len(allPPChrPos[chromosome])}\t{wdNoGene}\t{totWd}")
    
    # Total windows can be zero if the chromosome size is equal to window size. 
    # To avoid a division by 0 error, if total possible windows is 0, 
    # just assign 1 to the probability because all orthologs will be in that window anyway.
                    if totWd != 0:
                        pONEgenesWithinWindowChr = 1 - wdNoGene / totWd
                    else:
                        pONEgenesWithinWindowChr = 1
                   # print(pONEgenesWithinWindowChr)

                    
            pOFNoGeneInWindow = totalWdsWithoutAnyGenes / totalWdsonALLChromosomes
            #print(f"p of no gene in window = {pOFNoGeneInWindow}")

# Get the probability of at least this gene in at least one window (1 - p of no gene in window)
            pOfONEGeneInWindow = 1 - pOFNoGeneInWindow
            #print(f"p of at least ONE gene in window for = {pOfONEGeneInWindow}\n")
           
			# Calculate the summation of all the probabilities and 1-p to calculate geometric mean later on --- TO EXCLUDE ANCHOR I SUBTRACT THE log(P) OF ANCHOR LATER ON FROM THESE SUM
			# Since log (0) gives error ignore zeros. Prob. can be zero if there is just 1 gene on a chromosome that is smaller than window. e.g. chromosome- Y
			# Also prob. can be zero if there are orthologs on all chromosomes and not even a single window can be placed without any gene. Ignore in that case too
            import math

            if 0 < pOfONEGeneInWindow < 1:
                sum_log_p += math.log(pOfONEGeneInWindow)
                sum_log_1_minus_p += math.log(1 - pOfONEGeneInWindow)
            # Now I have to check which one of them exist in the block and push the probability at the end.
			# If there is no ortholog on any PP chr i.e. $_, i figured perl won't go in this loop. So only chromosomes having an ortholog are listed here.
			# If there are multiple polyploid orthologs on this gene, it doesn't matter as the probability will be same for them because they share same outgroup.
            
            for pChr in allPPChrPos.keys():
                if len(allPPChrPos[pChr]) == 0:
                    continue

                for ppPos in allPPChrPos[pChr]:
                
                    #print(f"*{og_anc_chr}\t{position}\t{pChr}\t{ppPos}")
                    
                    for blockLine in genes_in_this_block:
                
                        #print(blockLine)
                        
                        if blockLine.strip() == f"{og_anc_chr}\t{position}\t{pChr}\t{ppPos}":
                            #print(f"->{og_anc_chr}\t{position}\t{pChr}\t{ppPos}\t{pOfONEGeneInWindow}\t*{pONEgenesWithinWindowChr}*")
                            
                            genes_in_block_with_probability.append(f"{og_anc_chr}\t{position}\t{pChr}\t{ppPos}\t{pONEgenesWithinWindowChr}\t{pOfONEGeneInWindow}")
                            
                            og_genes_with_orth_in_this_wd[position] = pOfONEGeneInWindow  # To get unique outgroup genes having orthologs ON ENTIRE PP GENOME in this window - push the og position in key.
                            og_genes_with_orth_in_this_wd_chr[position] = pONEgenesWithinWindowChr  # To get unique outgroup genes having orthologs ON THIS CHROMOSOME OF PP GENOME in this window - push the og position in key.
                            
                            if f"{og_anc_chr}\t{position}\t{pChr}\t{ppPos}" in anchor:
                                p_anchor_genome = pOfONEGeneInWindow
                                p_anchor_chr = pONEgenesWithinWindowChr
                                #print(f"*{og_anc_chr}\t{position}\t{pChr}\t{ppPos}\t{pOfONEGeneInWindow}")    
#			print "\n----------------------------\n";
                                
                        #print("\n".join(genes_in_this_block))
        #print(uniqueGenesWithOrth)

#-- Calculate probability of finding all the unique Og genes in this window randomly ACROSS ENTIRE GENOME
    for key in og_genes_with_orth_in_this_wd.keys():
        # Here ignore the 0 values for individual genees. It can happen if there is just one ortholog on a chr/scaff which is smaller than the window size e.g. Y in human  # To ignore values that are 0
        #print(f"{key} => {og_genes_with_orth_in_this_wd[key]}")
        if og_genes_with_orth_in_this_wd[key] > 0:
            p_all_genes_within_window *= og_genes_with_orth_in_this_wd[key]
	#-- Calculate probability of finding all the unique Og genes in this window randomly ON THIS CHROMOSOME ONLY
    for key in og_genes_with_orth_in_this_wd_chr.keys():
        if og_genes_with_orth_in_this_wd_chr[key] > 0:
            #print(f"===> {key} => {og_genes_with_orth_in_this_wd_chr[key]}")
            p_all_genes_within_window_chr *= og_genes_with_orth_in_this_wd_chr[key]

    if p_anchor_chr is not None and p_anchor_chr > 0:   # Remove the probability for anchor
        p_all_genes_within_window_chr /= p_anchor_chr
    # To exclude P of anchor divide by anchor probability for joint probability where i multiply
	# And subtract of log values where I add.
	# Now this could be problematic as p for anchor can be zero if anchor is on chr smaller than PP wd
	# Also if there are so many genes that its not possible to find a window without any gene ******************* THINK ABOUT IT ALSO IF I NEED TO IGNORE THAT ???
	# So i only remove ancchor if its P is greater than 0. Else i dont remove it.
    
    import math
    if p_anchor_genome is not None and 0 < p_anchor_genome < 1:
        p_all_genes_within_window /= p_anchor_genome  # Since I have to EXCLUDE ANCHOR, divide by anchor probability for joint probability calculation
        sum_log_p -= math.log(p_anchor_genome)  # Deduct the log(P) for anchor from this sum 
        sum_log_1_minus_p -= math.log(1 - p_anchor_genome)  # Deduct the log(1 - P) for anchor from this sum
    #print(f"* P of all genes within window: {p_all_genes_within_window}")
    
	#my $pAllGenesWithinWindow = sprintf("%.3E", $pAllGenesWithinWindow);   # change to scientific notation

    #print(f"P of Anchor = {p_anchor_genome}")
    #print(f"SumLogP: {sum_log_p}, SumLog1MinusP: {sum_log_1_minus_p}")
    
    import math
    N = uniqueGenesWithOrth - 1  # Total outgroup genes having orthologs in PP genome-- EXCLUDING ANCHOR
    k = len(og_genes_with_orth_in_this_wd.keys()) - 1  # Total outgroup genes in this window having orthologs in PP window-- EXCLUDING ANCHOR
    LogpBar = sum_log_p / N  # Geometric mean of log P for entire window
    P1 = math.exp(LogpBar)
    LogOneMinusPBar = sum_log_1_minus_p / N  # Geometric mean of log 1-P for entire window
    P2 = 1 - math.exp(LogOneMinusPBar)
    #print(f"SumLogP: {sum_log_p}, SumLog1MinusP: {sum_log_1_minus_p}")
   #print(f"P1: {P1}, P2: {P2}")
    


    # Calculate Probability of observing any k number of genes within this window -- EXCLUDING ANCHOR	
    pAny_k_GenesWithinWindow_1 = None  # From log P
    pAny_k_GenesWithinWindow_2 = None  # From log 1-P

    # P from log approximation which Herv� told me to do earlier ---------- OLDER
	#$pAny_k_GenesWithinWindow_1 = (-1 * $k * (log($k/$N) - log($P1))) - (($N-$k) * (log(($N-$k)/$N) - log(1-$P1)));
	#$pAny_k_GenesWithinWindow_2 = (-1 * $k * (log($k/$N) - log($P2))) - (($N-$k) * (log(($N-$k)/$N) - log(1-$P2)));


	#--- Print probability in new anchor file, rest of the probabilities will be printed below when i evaluate the P for all possible K's
    with open('ANCPROB', 'a') as f:
        f.write(f"{anchor}\t{p_all_genes_within_window_chr}\t")
        f.write(f"{p_all_genes_within_window}\t")

#	print ANCPROB exp($pAny_k_GenesWithinWindow_1),"\t";
#	print ANCPROB exp($pAny_k_GenesWithinWindow_2),"\n";
	
	#--- Print probability of block in new block file, sameway rest of the probabilities will be printed below when i evaluate the P for all possible K's
    with open('PROB', 'a') as f:
        f.write(f">{anchor}\t{p_all_genes_within_window_chr}\t")
        f.write(f"{p_all_genes_within_window}\t")
#	print PROB exp($pAny_k_GenesWithinWindow_1),"\t";
#	print PROB exp($pAny_k_GenesWithinWindow_2),"\n";
	
#	foreach (@genesInBlockWithProbability){print PROB "$_\n";}
    # print(f"N = {N}\nk = {k}\n")
    # print(f"Sum log(p) = {sum_log_p}\n")
    # print(f"Log P bar = Sum_log(P)/N = {LogpBar}\n")
    # print(f"P1 = exp (Log P bar) = {P1}\n----------\n")
    # print(f"Sum log(1-p) = {sum_log_1_minus_p}\n")
    # print(f"Log 1 - Pbar = Sum log(1-p)/N = {LogOneMinusPBar}\n")
    # print(f"Exp (Log 1 - Pbar) = {math.exp(LogOneMinusPBar)}\n")
    # print(f"P2 = exp (1 - log(1 - Pbar)) = {P2}\n")
    
    sumkP1 = 0
    sumkP2 = 0  # Sum of P1 and P2 for all Ks

    P_MoreThanEqualK_1 = 0
    P_MoreThanEqualK_2 = 0  # P of finding any K or more orthologs randomly

    from math import comb, exp
    from decimal import Decimal, getcontext
    getcontext().prec = 100 
    
    for kLocal in range(k + 1):  # Iterate from 0 to k (inclusive)
    
        n = Decimal(N)  # Get the big float from N using Decimal for high precision
        # n = N  # Alternative, if no big float is needed

        N_C_k = comb(int(n), kLocal)  # Calculate "N choose k" using Python's comb function from the math module
        #print(f"{kLocal}\t{N}\t{N_C_k}")
   
        P1k = Decimal(P1) ** kLocal
        P2k = Decimal(P2) ** kLocal
        #print(f"\t{P1k}\t{P2k}")
        #print(f"P1: {P1}, P2: {P2}, N: {N}, kLocal: {kLocal}")
        
        OneMinusP1k = (Decimal(1) - Decimal(P1)) ** (Decimal(N) - kLocal)
        OneMinusP2k = (Decimal(1) - Decimal(P2)) ** (Decimal(N) - kLocal)
        
        #print(f"\t{OneMinusP1k}\t{OneMinusP2k}")
       
        # Convert to Decimal for high precision
        x1 = Decimal(N_C_k)
        x2 = Decimal(P1k)
        x3 = Decimal(OneMinusP1k)
        x4 = Decimal(P2k)
        x5 = Decimal(OneMinusP2k)

        # x1 = N_C_k
        # x2 = P1k
        # x3 = OneMinusP1k
        # x4 = P2k
        # x5 = OneMinusP2k

        # Set precision if needed (not directly applicable in Python without context)
        #print(f"\n->{x1}\n{x2}\n{x3}")
        
        testP1 = x1 * x2 * x3  # That's the P for the current k
        testP2 = x1 * x4 * x5

        #print(f"* {testP1}\n* {testP2}")
  
        # testP1 = f"{testP1:.3E}"
        # testP2 = f"{testP2:.3E}"
        # print(f"-> {testP1}\t{testP2}")

        # Sum P1 and P2's
        sumkP1 += testP1
        sumkP2 += testP2

        if kLocal == (k - 1):  # For all P's and sum of P below actual K
            
            P_MoreThanEqualK_1 = 1 - sumkP1
            P_MoreThanEqualK_2 = 1 - sumkP2
            #print(f"{kLocal}\n{sumkP1}\n{sumkP2}\n\n")
           
            with open('ANCPROB', 'a') as ancprob_file:
                ancprob_file.write(f"{P_MoreThanEqualK_1}\t{P_MoreThanEqualK_2}\t")  # Print P of finding any K or more genes in the block

            with open('PROB', 'a') as prob_file:
                prob_file.write(f"{P_MoreThanEqualK_1}\t{P_MoreThanEqualK_2}\t")  # Print P of finding any K or more genes in the block

        if k == kLocal:

            with open('ANCPROB', 'a') as ancprob_file:
                ancprob_file.write(f"{testP1}\t{testP2}\n")  # P of finding exactly any K gene randomly in the block

            with open('PROB', 'a') as prob_file:
                prob_file.write(f"{testP1}\t{testP2}\n")  # P of finding exactly any K gene randomly in the block
                for gene in genes_in_block_with_probability:
                    prob_file.write(f"{gene}\n")
            #print(f"-> {testP1}\t{testP2}\n")
        #print(f"{kLocal}\t{testP1}\t{testP2}\t{sumkP1}\t{sumkP2}\n")
        
#my $sumkP1_Formatted = sprintf("%.3E", $sumkP1);
	#my $sumkP2_Formatted = sprintf("%.3E", $sumkP2);
	
   
		

# Timing and output formatting
import time

startTime = time.time()

elapsed_time = time.time() - startTime
print("...[DONE] in ", end="")

if elapsed_time >= 3600:
    print(f"{elapsed_time / 3600:.2f} hours")
elif elapsed_time >= 60:
    print(f"{elapsed_time / 60:.2f} minutes")
else:
    print(f"{int(elapsed_time)} seconds")