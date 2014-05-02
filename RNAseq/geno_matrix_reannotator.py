start = open("genotypes.matrix","r")
end = open("genotypes.matrix.update","w")
for line in start:
    line = line.rstrip("\n")
    if "SampleID" in line:
        print >> end, line + "\tfindiv\tfc\tlane"
        continue
    findiv = line.split()[0].split(".")[0]
    fc = line.split()[0].split(".")[1][0:2]
    if fc[0] == 0:
        fc = fc[1]
    lane = line.split()[0].split(".")[1][3:4]
    print >> end, line + "\t" + str(findiv) + "\tfc" + str(fc) + "\tlane_" + str(lane)
start.close()
end.close()