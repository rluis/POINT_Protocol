import re
from typing import List

import pyensembl


class GenomicFeatures:
    """
    A set of functions to manipulate genomic features

    """

    @staticmethod
    def readBedFile(BED_PATH: str) -> list:
        """
        Read BED from file returning a nested list.
        :param BED_PATH: Bed file Path
        :return: Nested List of Bed Entries
        """

        # Read BED from file
        bedFile = open(BED_PATH)
        bedEntries = [row.strip().split() for row in bedFile.readlines() if not row.startswith("#")]
        bedFile.close()

        return bedEntries

    @staticmethod
    def writeBedFile(bedList: List[list], BED_PATH: str) -> None:
        """
        Write bedList to file.
        :param bedList: A list of lists of BED entries
        :param BED_PATH: Bed file path to write
        :return: None
        """

        # Read BED from file
        bedFile = open(BED_PATH, "w")

        # Write BED file
        bedFile.write('\n'.join('\t'.join(map(str, row)) for row in bedList))

        # Close file
        bedFile.close()


    @staticmethod
    def sortBedList(bedList: List[list]) -> List[list]:
        """
        Sorts a list of BED lists.
        Sorts by: chr1, chr2, chr3... chr10,chr11,chr12... chr20, chr21

        :param bedList: A list of lists of BED entries
        :return: A sorted bedList
        """

        def keySortFuncBedFile(s):
            return [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s[0])], int(s[1])

        return sorted(bedList, key=keySortFuncBedFile)

    @staticmethod
    def createGeneBED(ensemblDB: pyensembl, geneID: str) -> list:
        """
        Creates a gene bed entry.

        :param ensemblDB: a pyensembl object.
        :param geneID: the gene ID for which the function will find the coordinates.
        :return: a BED row (in a list format).
        """
        chr = ensemblDB.gene_by_id(geneID).contig
        start = str(ensemblDB.gene_by_id(geneID).start)
        end = str(ensemblDB.gene_by_id(geneID).end)
        name = geneID
        score = "0"
        strand = str(ensemblDB.gene_by_id(geneID).strand)
        return [chr,start,end,name,score,strand]

    @staticmethod
    def createTransBED(ensemblDB: pyensembl, transcriptID: str) -> list:
        """
        Creates a transcript bed entry.

        :param ensemblDB: a pyensembl object.
        :param transcriptID: the transcript ID for which the function will find the coordinates.
        :return: a BED row (in a list format).
        """
        chr = ensemblDB.transcript_by_id(transcriptID).contig
        start = str(ensemblDB.transcript_by_id(transcriptID).start)
        end = str(ensemblDB.transcript_by_id(transcriptID).end)
        name = transcriptID
        score = "0"
        strand = str(ensemblDB.transcript_by_id(transcriptID).strand)
        return [chr,start,end,name,score,strand]

    @staticmethod
    def transcriptID_exonsIDs(ensemblDB: pyensembl, transcriptID: str) -> list:
        """
        Returns the a list of all exonIDs that are associated with the given TranscriptID.

        :param ensemblDB: a pyensembl object.
        :param transcriptID: the transcriptID to access
        :return: Returns a list of all exonIDs associated
        """
        return ensemblDB.exon_ids_of_transcript_id(transcriptID)

    @staticmethod
    def createExonBED(ensemblDB: pyensembl, exonID: str, transcriptID: str) -> list:
        """
        Creates an exon bed entry.

        :param ensemblDB: a pyensembl object.
        :param exonID: the exon ID for which the function will find the coordinates.
        :param transcriptID: the transcript ID for which the function will find the coordinates.
        :return: a BED row (in a list format).
        """
        chr = ensemblDB.exon_by_id(exonID).contig
        start = str(ensemblDB.exon_by_id(exonID).start)
        end = str(ensemblDB.exon_by_id(exonID).end)
        name = "{}_{}".format(exonID, transcriptID)
        score = 0
        strand = str(ensemblDB.exon_by_id(exonID).strand)
        return [chr, start, end, name, score, strand]

    @staticmethod
    def transcriptID_exonsBedList(ensemblDB: pyensembl, transcriptID: str) -> List[list]:
        """
        Returns a BedList of all exons of a given transcriptID.
        Feature name is the ExonID. addTransID and addGeneID can be used to increment the respective IDs.

        :param ensemblDB: a pyensembl object.
        :param transcriptID: the transcriptID to access
        :return: Retruns a bedList of all exons associated to a transcriptID.
        """
        tmpList = []
        for exonID in GenomicFeatures.transcriptID_exonsIDs(ensemblDB=ensemblDB, transcriptID=transcriptID):

            tmpList.append(GenomicFeatures.createExonBED(ensemblDB=ensemblDB,
                                                        exonID=exonID,
                                                        transcriptID=transcriptID))

        return tmpList

    @staticmethod
    def createIntronBedEntry(ensemblDB: pyensembl, exonID1: str, exonID2: str, transcriptID: str) -> list:
        """
        Creates a intron bedEntry from two exonIDs and transcriptID.

        :param ensemblDB: An object created of ensemblDB
        :param exonID1: one exonID to access. Not necessary to give in a particular order. Exons are sorted inside by coordinate.
        :param exonID2: second exonID to access. Not necessary to give in a particular order. Exons are sorted inside by coordinate.
        :param transcriptID: TranscriptID of the desired exon.
        :return: A bedEntry of a intron between the given 2 exonIDs of the transcriptID.
        """

        # Checks if given exons are in the same transcripts
        exonsOfTheTransID = GenomicFeatures.transcriptID_exonsIDs(ensemblDB, transcriptID)
        if exonID1 not in exonsOfTheTransID or exonID2 not in exonsOfTheTransID:
            raise ValueError("At least one of the given exons are not associated with the trasncriptID")

        # Crete bedEntries of the exons
        exonBedEntry_1 = GenomicFeatures.createExonBED(ensemblDB, exonID1, transcriptID)
        exonBedEntry_2 = GenomicFeatures.createExonBED(ensemblDB, exonID2, transcriptID)
        exonBedEntry_1_sort, exonBedEntry_2_sort = GenomicFeatures.sortBedList([exonBedEntry_1, exonBedEntry_2])

        # Build intron name
        featureName = "{}_{}_{}".format(exonBedEntry_1_sort[3].split("_")[0], exonBedEntry_2_sort[3].split("_")[0], transcriptID)

        # Create BedEntry Intron
        return [exonBedEntry_1_sort[0],
                str(exonBedEntry_1_sort[2]),
                str(exonBedEntry_2_sort[1]),
                featureName,
                "0",
                exonBedEntry_1_sort[5]]

    @staticmethod
    def transcriptID_intronsBedList(ensemblDB: pyensembl, transcriptID: str) -> List[list]:
        """
        Returns a BedList of all introns of a given transcriptID.
        Feature name is the firstExonID_secondExonID, sorted by coordinates.
        addTransID and addGeneID can be used to increment the respective IDs.

        :param ensemblDB: An object created of ensemblDB
        :param transcriptID: the transcriptID to access
        :param addTransID: if True, adds the transID to the bedEntry name (ENSE0000_ENSE0001_ENST0000)
        :param addGeneID: if True, adds the GeneID to the bedEntry name (ENSE0000_ENSE0001_ENSG0000,
                            or ENSE0000_ENSE0001_ENST0000_ENSG0000 if combined with addTransID)
        :return: Retuns a list of introns associated a transcriptID
        """

        # Get all exons of the transID
        exonsBedList = GenomicFeatures.transcriptID_exonsBedList(ensemblDB=ensemblDB,transcriptID=transcriptID)
        exonsBedListSorted = GenomicFeatures.sortBedList(exonsBedList)

        # If intronless returns a empty list
        if len(exonsBedListSorted) <= 1:
            return []

        # If not intronless builds intron features
        tmpList = []
        for enum, exonEntry in enumerate(exonsBedListSorted):
            if enum == 0:
                prevEntry = exonEntry
            else:
                tmpList.append(GenomicFeatures.createIntronBedEntry(ensemblDB, prevEntry[3].split("_")[0],
                                                                 exonEntry[3].split("_")[0], transcriptID=transcriptID))
                prevEntry = exonEntry

        return tmpList
