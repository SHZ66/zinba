import CountMap

print "track type=wiggle_0 name=\"HG18 Alignability (36bp)\" desc=\"HG18 Alignability (36bp)\" visibility=full"

for i in xrange(1,26):
        if i <= 22:
                chr = str(i)
        elif i == 23:
                chr = "M"
        elif i == 24:
               chr = "X"
        elif i == 25:
                chr = "Y"

        nicechr = "chr" + chr
        filename = "chr" + chr + "b.out"
        cm=CountMap.CountMap(filename)
        print "fixedStep chrom=" + nicechr + " start=1 step=1"

        for j in xrange(1,21000000):        
                flag = 0
                try:
                        x = cm.cnt(j)			
                except ValueError:
                        flag = 1
                if flag == 0:
                        print str(x)

