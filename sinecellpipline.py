# Autor: chentao
# E-mail: hmcctt@163.com
# Created on: 2019-07-18

import subprocess
import sys
import os
import re
import argparse

parse = argparse.ArgumentParser(description="Sing cell sequencing project.",
                                epilog="Created on 2019-07-08.")
parse.add_argument('--level', '-l', choices=['expression', 'variant'],
                   required=True, help='choice your pipline with variant or expression')
parse.add_argument('--pipline', '-p', choices=['star', 'rsem'], default='star',
                   help='choose pipline to mapping reads')
parse.add_argument('--thread', '-t', type=int, default=6,
                   help='threads')
parse.add_argument('--input', '-i', required=True, help='your input fires')
parse.add_argument('--output', '-o', required=True, help='your out put directory')

args = parse.parse_args()


class fastqdump():
    
    def __init__(self, input='/input/testing', output='/output/testing',
                 thread = 6):
        self.input = input
        self.ouput = output
        self.thread = thread
        self.filename = os.path.basename(input).replace('.sra', '')

    def runstatus(self, command):
        a = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if a.returncode != 0:
            print(a.stderr)
            sys.exit(120)

    def RunFastqdump(self):
        print('RunFastqdump')
        commands = ['fastq-dump', '--split-3', '-O']
        commands.append(self.ouput)
        commands.append(self.input)
        self.runstatus(commands)

    def Trimmomatic(self, left):
        print('Trimmomatic')
        commands = ['java', '-jar', '/home/tchen/Chentao/software/Trimmomatic-0.38/trimmomatic-0.38.jar',
                    'PE', '-threads', str(self.thread)]
        commands.append(self.ouput + '/' + self.filename + '_1.fastq')
        commands.append(self.ouput + '/' + self.filename + '_2.fastq')
        commands.append(self.ouput + '/' + self.filename + '_1.clean.fastq')
        commands.append(self.ouput + '/' + self.filename + '_1.unparied.fastq')
        commands.append(self.ouput + '/' + self.filename + '_2.clean.fastq')
        commands.append(self.ouput + '/' + self.filename + '_2.unparied.fastq')
        commands.append('HEADCROP:' + str(left))
        self.runstatus(commands)

    def Star(self, genome, type = 'expression'):
        print('Star')
        commands = ['STAR', '--outFilterType BySJout', '--outFilterMultimapNmax 20',
                    '--alignSJoverhangMin 8', '--alignSJDBoverhangMin 1',
                    '--outFilterMismatchNmax 999', '--outFilterMismatchNoverReadLmax 0.04',
                    '--alignIntronMin 20', '--alignIntronMax 1000000',
                    '--alignMatesGapMax 1000000', '--outSAMstrandField intronMotif',
                    '--outSAMtype BAM SortedByCoordinate'
                    ]
        commands.append('--outSAMattrRGline ' + "ID:" + self.filename + " SM:" + self.filename +
                        ' PL:illumina' + ' LB:' + self.filename + ' PU:' + self.filename)
        commands.append('--runThreadN ' + str(self.thread))
        commands.append('--readFilesIn ' +
                        self.ouput + '/' + self.filename + '_1.clean.fastq' + " " +
                        self.ouput + '/' + self.filename + '_2.clean.fastq')
        commands.append('--genomeDir ' + genome)
        if type != 'expression':
            commands.append('--twopassMode Basic')
        commands.append('--outFileNamePrefix ' + self.ouput + '/' + self.filename + '.')
        self.runstatus(commands)

    def Picard(self):
        print('Picard')
        commands = ['java', '-jar',
                    '/home/tchen/Chentao/picard.jar', 'MarkDuplicates',
                    'REMOVE_DUPLICATES=true']
        commands.append('I=' + self.ouput + '/' + self.filename + '.Aligned.sortedByCoord.out.bam')
        commands.append('O=' + self.ouput + '/' + self.filename + '.md.bam')
        commands.append('M=' + self.ouput + '/' + self.filename + '.md.metrics.txt')
        commands.append('CREATE_INDEX=true')
        self.runstatus(commands)

    def HtseqCount(self):
        print('htseq-count')
        commands = ['htseq-count', '-m', 'intersection-nonempty',
                    '-s', 'no', '-i', 'gene_id', '-f', 'bam']
        commands.append(self.ouput + '/' + self.filename + '.md.bam')
        commands.append('/home/tchen/Chentao/gliomas_scRNA/genome' +
                        '/gencode.v19.annotation.gtf')
        outputfiles = self.ouput + '/' + 'readcount' + '/' + self.filename + '.txt'
        a = subprocess.run(commands, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if a.returncode != 0:
            print(a.stderr)
            sys.exit(120)
        else:
            print(a.stdout.decode('utf8'), file=open(outputfiles, 'w'))

    def Rsem(self):
        print("Rsem")
        commands = ['/home/tchen/Chentao/software/RSEM-master/rsem-calculate-expression',
                    '--paired-end']
        commands.extend(['-p', str(self.thread)])
        commands.append(self.ouput + '/' + self.filename + '_1.fastq')
        commands.append(self.ouput + '/' + self.filename + '_2.fastq')
        commands.append('/home/tchen/Chentao/gliomas_scRNA/genome/bowtie/human')
        commands.append(self.ouput + '/' + self.filename)
        self.runstatus(commands)
        a = subprocess.run(['mv', self.ouput + '/' + self.filename + '.genes.results',
                            self.ouput + '/readcount/'],
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE)
        if a.returncode != 0:
            print(a.stderr)
            sys.exit(120)

if os.path.exists(args.output):
    if os.path.exists(args.output + '/' + 'readcount'):
        pass
    else:
        os.mkdir(args.output + '/' + 'readcount')
else:
    print('output dir is not exits.')

testing = fastqdump(args.input, args.output, args.thread)
testing.RunFastqdump()


fastfilename = os.path.basename(args.input).replace('.sra', '')
if args.level == 'expression':
    if args.pipline == 'star':
        files = open(args.output + '/' + fastfilename + '_1.fastq', 'r')
        seq_len = re.search(r'length=([0-9]+)', files.readline()).group(1)
        files.close()

        if int(seq_len) > 90:
            left = 20
            genome = '/home/tchen/Chentao/gliomas_scRNA/genome/star/human_l80'
        else:
            left = 5
            genome ='/home/tchen/Chentao/gliomas_scRNA/genome/star/human_l60'
        testing.Trimmomatic(left)
        testing.Star(genome)
        testing.Picard()
        testing.HtseqCount()
    else:
        testing.Rsem()
