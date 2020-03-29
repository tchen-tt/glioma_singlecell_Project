# Created on: 2019-12-09
# E-mail: hmcctt@163.com
import os
import argparse
import subprocess

parse = argparse.ArgumentParser(description="Caculate the gene exression level",
                                epilog="Created on 2019-07-08.")
parse.add_argument('--thread', '-t', type=int, default=6,
                   help='threads')
parse.add_argument('--input', '-i', required=True, help='your input fires')
parse.add_argument('--output', '-o', required=True, help='your out put directory')

args = parse.parse_args()



class RnaExpress():
    def __init__(self, input, output, threads = 6):
        self.input = input
        self.output = output
        self.threads = threads
        self.filename = os.path.basename(input).replace(".sra", "")

    def runstatus(self, command):
        a = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if a.returncode != 0:
            print(a.stderr)
            sys.exit(120)

    def RunFastqdump(self):
        print('RunFastqdump')
        commands = ['fastq-dump', '--split-3', '-O']
        commands.append(self.output)
        commands.append(self.input)
        self.runstatus(commands)

    def Trimmomatic(self, left):
        print('Trimmomatic')
        commands = ['java', '-jar', '/home/tchen/Chentao/software/Trimmomatic-0.38/trimmomatic-0.38.jar',
                    'SE', '-threads', str(self.threads)]
        commands.append(self.output + '/' + self.filename + '.fastq')
        commands.append(self.output + '/' + self.filename + '.clean.fastq')
        commands.append('HEADCROP:' + str(left))
        self.runstatus(commands)


    def Rsem(self):
        print("Rsem")
        commands = ['/home/tchen/Chentao/software/RSEM-master/rsem-calculate-expression',
                    '--bowtie2']
        commands.extend(['-p', str(self.threads)])
        commands.append(self.output + '/' + self.filename + '.clean.fastq')
        commands.append('/home/tchen/Chentao/gliomas_scRNA/genome/bowtie2/human/GRch37.p13.genecode')
        commands.append(self.output + '/' + self.filename)
        self.runstatus(commands)
        a = subprocess.run(['mv', self.output + '/' + self.filename + '.genes.results',
                            self.output + '/readcount/'],
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE)
        if a.returncode != 0:
            print(a.stderr)
            sys.exit(120)

if __name__ == "__main__":
    if os.path.exists(args.output + '/' + 'readcount'):
        pass
    else:
        os.mkdir(args.output + '/' + 'readcount')

    start = RnaExpress(args.input, args.output, args.thread)
    start.RunFastqdump()
    start.Trimmomatic(15)
    start.Rsem()

