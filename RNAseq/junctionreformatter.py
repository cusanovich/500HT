"""Take in a sam file and write out a reformatted file."""
import sys

f = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')
for line in f:
    if line == "":
        continue
    else:
        line = line.strip().split('\t')
        if line[1] != '0' and line[1] != '16':
            continue
        elif int(line[3]) > 50:
            continue
        elif int(line[3]) < 2:
            continue
        else:
            liner = ['','','','','','']
#            print line
#            break
            liner[0] = line[2].split('_')[0]
            liner[1] = str(int(line[2].split('_')[1]) - 51 + int(line[3]))
            liner[2] = line[2].split('_')[1]
            liner[3] = line[0]
            liner[4] = line[4]
            if line[1] == '0':
                liner[5] = '+'
            elif line[1] == '16':
                liner[5] = '-'
#            print liner
#            break
            print >> outfile, '\t'.join(liner)
f.close()
outfile.close()