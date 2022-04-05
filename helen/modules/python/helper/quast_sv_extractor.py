import argparse
"""
HELPER SCRIPT TO COUNT QUAST MISASSEMBLIES THAT HAPPEN OUTSIDE KNOWN SVs.
"""

autosomes = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY'] + ['chrX_hg002', 'chrY_hg002'] + ["GWHBDNP00000001.1_OriSeqID_chr1_Len_32659241","GWHBDNP00000002.1_OriSeqID_chr2_Len_22560461","GWHBDNP00000003.1_OriSeqID_chr3_Len_26161332","GWHBDNP00000004.1_OriSeqID_chr4_Len_22250686","GWHBDNP00000005.1_OriSeqID_chr5_Len_30093473"]

def read_quast_file(file, tigList):
    """
    Extracts misassembly regions from QUAST file, usually located in /QUAST_output/contig_reports/all_alignments_$$.tsv
    :param file:
    :return:
    """
    total_count = 0
    relocation_count = 0
    translocation_count = 0
    inversion_count = 0
    misassemblies = []
    isError = False
    with open(file) as f:
        prev_line = ''
        for line in f:
            line = line.rstrip().replace(',', '')
            splits = line.split('\t')
            if isError == True:
                #print("Ready to add error on the other side and my current ref start and end are %s and %s and on the tig its %s and %s"%(s_ref, e_ref, s_con, e_con)) 
                if (int(s_con) < int(e_con)):
                   s=e_ref
                else:
                   s=s_ref
 
                s_ref, e_ref, s_con, e_con, ref, con, idn, ambi, bg = line.split('\t')
                isError = False
                #print("Aadding error after type in %s from %s to %s"%(ref, s_ref, e_ref))
                if errorType == 'relocation':
                   # now we can add this because it's giving us the other coordinate 
                   if len(tigList) == 0 or con in tigList:
                      if (int(s_con) < int(e_con)):
                         misassemblies.append([ref, min(s, s_ref), max(s, s_ref), errorType])
                         print("Adding relocation in fwd strand ", ref, min(s, s_ref), max(s, s_ref), splits[0])
                      else:
                         misassemblies.append([ref, min(s, e_ref), max(s, e_ref), errorType])
                         print("Adding relocation in rev strand ", ref, min(s, e_ref), max(s, e_ref), splits[0])
                      relocation_count += 1
                      total_count += 1
                elif errorType == 'translocation':
                   # we add both ends of a translocation and count it as two errors cause it has two locations on the reference right
                   if len(tigList) == 0 or con in tigList:
                      #print("Reading error for tig %s\n", con)
                      print("Adding second translocation ", ref, s_ref, e_ref, splits[0])
                      misassemblies.append([ref, s_ref, e_ref, errorType])
                      translocation_count += 1
                      total_count += 1
                elif errorType == 'inversion':
                   if len(tigList) == 0 or con in tigList:
                      #print("Reading error for tig %s\n", con)
                      print("Adding inversion ", ref, s_ref, e_ref, splits[0])
                      misassemblies.append([ref, s_ref, e_ref, errorType])
                      inversion_count += 1
                      total_count += 1
            errorType = splits[0].split(' ')[0]

            if splits[0].split(' ')[0] == 'relocation':
                isError=True
                s_ref, e_ref, s_con, e_con, ref, con, idn, ambi, bg = prev_line.split('\t')
                #print(ref, s_ref, e_ref, splits[0])
            elif splits[0] == 'translocation':
                isError=True
                s_ref, e_ref, s_con, e_con, ref, con, idn, ambi, bg = prev_line.split('\t')
                print("Adding first translocation ", ref, s_ref, e_ref, splits[0])
                if len(tigList) == 0 or con in tigList:
                   #print("Reading error for tig %s\n", con)
                   misassemblies.append([ref, s_ref, e_ref, splits[0]])
                   translocation_count += 1
                   total_count += 1
            elif splits[0] == 'inversion':
                isError=True
                s_ref, e_ref, s_con, e_con, ref, con, idn, ambi, bg = prev_line.split('\t')
                #print(ref, s_ref, e_ref, splits[0])
                # we will add the inversion later since that is reporting the actual inverted sequence on the reference
            prev_line = line

    print("#####--TOTAL MISSASSEMBLIES: REPORTED BY QUAST--#####")
    print("Total Misassemblies:\t", total_count)
    print("Total relocations:\t", relocation_count)
    print("Total translocations:\t", translocation_count)
    print("Total inversions:\t", inversion_count)
    print("#####################################################\n")
    # returns a list of misassemblies as a list of list where each element is [chr_name, st, end, type_of_missassmbly]
    return misassemblies

def read_tig_file(file):
    known_tigs = []
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            known_tigs.append(line.split()[0])
            # print(line.split('\t')[0:3])
    # returns a list of known svs as a list where each element is [chr_name, st, end]
    return known_tigs

def read_bed_file(file):
    known_svs = []
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            known_svs.append(line.split('\t')[0:3])
            # print(line.split('\t')[0:3])
    # returns a list of known svs as a list where each element is [chr_name, st, end]
    return known_svs


def count_miassemblies_in_autosomes(misassemblies):
    total_count = 0
    relocation_count = 0
    translocation_count = 0
    inversion_count = 0

    for chr, st, end, ms_type in misassemblies:
        if chr in autosomes:
            if ms_type == 'relocation':
                relocation_count += 1
                total_count += 1
            elif ms_type == 'translocation':
                translocation_count += 1
                total_count += 1
            elif ms_type == 'inversion':
                inversion_count += 1
                total_count += 1
        else:
            print(chr, st, end, ms_type)

    print("#####--TOTAL MISSASSEMBLIES: IN AUTOSOMES--#####")
    print("Total Misassemblies:\t", total_count)
    print("Total relocations:\t", relocation_count)
    print("Total translocations:\t", translocation_count)
    print("Total inversions:\t", inversion_count)
    print("################################################\n")


def count_misassemblies_not_overlapping_with_svs(known_svs, misassemblies):
    total_count = 0
    total_base_count = 0
    relocation_count = 0
    translocation_count = 0
    inversion_count = 0

    for chr_ms, st_ms, end_ms, ms_type in misassemblies:
        if chr_ms not in autosomes:
            continue

        overlaps = False
        for chr_sv, st_sv, end_sv in known_svs:
            #print("Comparing sv %s %s %s to  known %s %s %s and booleans are %s %s %s"%(chr_ms, st_ms, end_ms, chr_sv, st_sv,  end_sv, (chr_ms == chr_sv), (int(st_ms) <= int(end_sv)), (int(st_sv) < int(end_ms))))
            # x1 <= y2 && y1 <= x2
            if chr_ms == chr_sv and int(st_ms) <= int(end_sv) and int(st_sv) <= int(end_ms):
                #print(chr_ms, st_ms, end_ms, ms_type, 'overlaps with sv', chr_sv, st_sv, end_sv)
                overlaps = True
                break
            if overlaps:
                break

        if overlaps is False:
            print(chr_ms, st_ms, end_ms, ms_type)
            total_count += 1
            total_base_count += (int(end_ms) - int(st_ms) + 1)
            if ms_type == 'relocation':
                relocation_count += 1
            elif ms_type == 'translocation':
                translocation_count += 1
            elif ms_type == 'inversion':
                inversion_count += 1

    print("#####--MISSASSEMBLIES OUTSIDE SVs--#####")
    print("Total Misassemblies:\t", total_count)
    print("Total relocations:\t", relocation_count)
    print("Total translocations:\t", translocation_count)
    print("Total inversions:\t", inversion_count)
    print("Total Bases:\t", total_base_count/ 1000000000)
    print("########################################")


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="Reports QUAST misassemblies that happen outside the known SVs.")
    parser.add_argument(
        "-q",
        "--quast_file",
        type=str,
        required=True,
        help="Path to the QUAST's TSV input file located in /QUAST_output/contig_reports/all_alignments_$$.tsv"
    )
    parser.add_argument(
        "-s",
        "--sv_file",
        type=str,
        required=True,
        help="Path to a bed file containing all known SVs for the sample."
    )
    parser.add_argument(
        "-c",
        "--cen_file",
        type=str,
        required=True,
        help="Path to a bed file containing all known centromeric region."
    )
    parser.add_argument(
        "-d",
        "--segdup_file",
        type=str,
        required=True,
        help="Path to a bed file containing all known segdup for the sample."
    )
    parser.add_argument(
        "-t",
        "--tig_file",
        type=str,
        required=True,
        help="Path to a list of tigs to include."
    )
    FLAGS, unparsed = parser.parse_known_args()
    # get regions from files
    tigsList = read_tig_file(FLAGS.tig_file)
    svs = read_bed_file(FLAGS.sv_file)
    centromeres = read_bed_file(FLAGS.cen_file)
    segdups = read_bed_file(FLAGS.segdup_file)

    centromeres_only = centromeres
    centromeres_and_segdups = centromeres + segdups
    m_assemblies = read_quast_file(FLAGS.quast_file, tigsList)

    # count stuff
    # count_miassemblies_in_autosomes(m_assemblies)
    print("OUTSIDE CENTROMERES")
    count_misassemblies_not_overlapping_with_svs(centromeres_only, m_assemblies)
    print("\nOUTSIDE CENTROMERES AND SEG DUPS")
    count_misassemblies_not_overlapping_with_svs(centromeres_and_segdups, m_assemblies)
    print("\nOUTSIDE CENTROMERES AND SV")
    count_misassemblies_not_overlapping_with_svs(centromeres + svs, m_assemblies)
    print ("\nOUTSIDE CENTROMERES AND SEG DUPS AND SVS")
    count_misassemblies_not_overlapping_with_svs(centromeres_and_segdups + svs, m_assemblies)
