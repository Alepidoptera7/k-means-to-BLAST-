import sys
import random
import math
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''

    def __init__(self, fname= 'read1.fa'):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            #header = ''
            #sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
                header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class k_means:
    """
    This class develops genomic categorization.
    """

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        self.read_dict = {}
        self.k = 9

        self.multi_read_cluster_holder_dict = {}
        self.centroids_dict = {}
        self.multi_read_centroid_change_list = [i for i in range(self.k)]
        self.seq_dict = {}
        self.sequence_assembly_dict = {}

        self.x_coords_dict = {}
        self.y_coords_dict = {}

        self.random_string = ''

    def rand_seq(self, seq):
        """Develops a random string the length of the read sequence for testing and comparison.

        Input: original sequence
        Output: random sequence of the same length
        """

        return ''.join(random.choices(['A','G','T','C'], k=len(seq)))


    def pr_matrix(self, head, seq):
        """

        Develops a dictionary in which each kmer key is associated with a list of unique (x,y) = (position, frequency) coordinate tuples.

        These values are attached to the given kmer and used at y-coordinates for Euclidian distance calculations.

        1) The sequence is parsed into 4-mers.
        2) The 4-mers are used for keys in the dictionary.
        3) The number of times each 4-mer occurs in the sequence is counted and divided by the number of possible 4mers, 256.


        Input: FASTA sequence
        Output: Dictionary with frequency Pr values attached to kmers.
        """

        #seq = self.rand_seq(seq)

        self.kmers = [seq[i:i + 4] for i in range(0, len(seq), 4)]
        self.seq_dict[head] = seq
        self.sequence_assembly_dict[head] = {}

        kmer_xy_coordinate_list = []
        for i in range(0, len(self.kmers)):
            if len(self.kmers[i]) == 4:
                kmer_xy_coordinate_list.append((i*4, self.kmers.count(self.kmers[i]) / 256))

        self.read_dict[head] = kmer_xy_coordinate_list

    def random_centroids(self):
        """Randomly selects k number of initial centroids by position in FASTA sequence.

        The initial centroids are attached to the read in the overarching read_dict.

        These centroids will be updated and used for cluster identification.

        Input: none
        Output: dictionary random centroids

        """

        for read in self.read_dict.keys():
            random_centroids = []

            for i in range(0, self.k):
                random_centroids.append(self.read_dict[read][random.randint(0, self.k)])

            self.centroids_dict[read] = random_centroids
            self.multi_read_cluster_holder_dict[read] = {i: [] for i in range(0, self.k)}

    def distance_calculation(self):
        """Euclidian distance is calculated between each centroid and all non-self points.

        Input: a list of points and centroids.
        Output: a list of distances between each centroid and all other points.
        """

        for read in self.read_dict.keys():
            temp_placement_dict = {i: [] for i in range(0, self.k)}

            for point in self.read_dict[read]:
                distance_list = []
                x2,y2 = point[0], point[1]

                for centroid in self.centroids_dict[read]:
                    x1, y1 = centroid[0], centroid[1]

                    #distance calculation between centroid and point
                    distance = math.sqrt(abs(x2 - x1) + abs(y2 - y1))
                    distance_list.append(distance)

                temp_placement_dict[distance_list.index(min(distance_list))].append(point)

            for i in range(self.k):
                self.multi_read_cluster_holder_dict[read][i] = temp_placement_dict[i]


    def centroid_recalculation(self):
        """Recalculates centroids by cluster averaging.

        Input: Current cluster information
        Output: List of adjusted centroids as component to self.read_dict[read][1]
        """

        for read in self.multi_read_cluster_holder_dict.keys():
            temp_centroid_list = []
            centroid_change_list = []
            for i in range(self.k):
                if len(self.multi_read_cluster_holder_dict[read][i]) > 0:
                    new_centroid_x = sum([tup[0] for tup in self.multi_read_cluster_holder_dict[read][i]]) / len(self.multi_read_cluster_holder_dict[read][i])
                    new_centroid_y = sum([tup[1] for tup in self.multi_read_cluster_holder_dict[read][i]]) / len(self.multi_read_cluster_holder_dict[read][i])

                    temp_centroid_list.append((new_centroid_x, new_centroid_y))
                    centroid_change_list.append(abs(new_centroid_x - self.centroids_dict[read][i][0]) + abs(new_centroid_y - self.centroids_dict[read][i][1]))

            self.multi_read_centroid_change_list = centroid_change_list
            self.centroids_dict[read] = temp_centroid_list


    def sequence_assembly(self):
        """Assembling the sequences represented by the predicted clusters.

        Input: Finished Clusters
        Output: Separated genomic sequences
        """

        for read in self.multi_read_cluster_holder_dict.keys():
            print(read)

            self.x_coords_dict[read] = {i: [] for i in range(0, self.k)}
            self.y_coords_dict[read] = {i: [] for i in range(0, self.k)}

            for i in self.multi_read_cluster_holder_dict[read].keys():
                cluster_sequence = ""
                self.sequence_assembly_dict[read][i] = ""

                x_coords_list = []
                y_coords_list = []

                for point in self.multi_read_cluster_holder_dict[read][i]:
                    cluster_sequence += self.seq_dict[read][point[0]: point[0] + 4]

                    x_coords_list.append(point[0])
                    y_coords_list.append(point[1])

                self.sequence_assembly_dict[read][i] = cluster_sequence

                self.x_coords_dict[read][i] = x_coords_list
                self.y_coords_dict[read][i] = y_coords_list

                if cluster_sequence != "":

                    print("cluster: ", i, " cluster size: ", "(", len(self.multi_read_cluster_holder_dict[read][i]), "/", len(self.read_dict[read]), ")", " portion of read: ",
                         len(self.multi_read_cluster_holder_dict[read][i]) / len(self.read_dict[read]))

                    print(self.sequence_assembly_dict[read][i])

                    self.blast(cluster_sequence)

                    print("_____________________________________________________________________________\n")
                    print("_____________________________________________________________________________\n")


    def cluster_graphing(self):
        """ The purpose of this function is to develop graphs of each computational cluster.
        Each subcluster will be represented in a different color, with all centroids marked similarly.

        Input: Dictionaries of x and y coordinates per cluster per read, dictionaries of centroids
        Output: Graphed clusters, each cluster in a distinct color in respect to each read.
        """

        for read in self.multi_read_cluster_holder_dict.keys():

            for i in self.multi_read_cluster_holder_dict[read].keys():

                xpoints_for_plotting = np.array(self.x_coords_dict[read][i])
                ypoints_for_plotting = np.array(self.y_coords_dict[read][i])

                color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:cyan', 'tab:purple', 'tab:brown',
                              'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:cyan', 'tab:purple', 'tab:brown']

                plt.scatter(xpoints_for_plotting, ypoints_for_plotting, color=color_list[i], marker='o')

                centroids_x_list = [i for i, j in self.centroids_dict[read]]
                centroids_y_list = [j for i, j in self.centroids_dict[read]]
                centroid_x_points_for_plotting = np.array(centroids_x_list)
                centroid_y_points_for_plotting = np.array(centroids_y_list)

                plt.scatter(centroid_x_points_for_plotting, centroid_y_points_for_plotting, color='k', marker='*')

            plt.title(label=read)
            plt.xlabel('kmer position')
            plt.ylabel('kmer frequency')

            plt.savefig(f"{read}.png")

        plt.show()


    def blast(self, cluster_sequence):
        """Develops a BLAST request for each clustered sequence via BioPython.

        Design taken from Biopython Manual.
        """
        E_VALUE_THRESH = 0.0000000000000000001

        result_handle = NCBIWWW.qblast("blastn", "nt", sequence=cluster_sequence, expect=1)

        blast_record = NCBIXML.read(result_handle)

        count = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    count += 1
                    print(f"sequence: {alignment.title}, length: {alignment.length}nt, e-value: {hsp.expect}")
                    #print(hsp.query[0:75] + "...")
                    #print(hsp.match[0:75] + "...")
                    #print(hsp.sbjct[0:75] + "...")

        print(f"There are {count} similar sequences in Blast output\n")


    def driver(self):
        """Drivers iterative processes"""

        self.random_centroids()

        self.distance_calculation()
        iteration = 0

        while sum(self.multi_read_centroid_change_list) > 0:
            self.centroid_recalculation()
            self.distance_calculation()

            iteration += 1

        print("total iterations of k-means clustering: ", iteration)

        self.sequence_assembly()
        self.cluster_graphing()


def main():

    class_access = k_means()

    myReader = FastAreader()

    for head, seq in myReader.readFasta():
        class_access.pr_matrix(head, seq)

    class_access.driver()

if __name__ == '__main__':
    main()