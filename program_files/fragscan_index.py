import sys

f = open(sys.argv[1])
fraggene = f.read().strip().split("\n")
f.close()

filetype = {
    "bed":{"Subject ID":3,"Percentage of identical matches":4,"Start of alignment in query":1,"End of alignment in query":2},
    "out":{"Subject ID":1,"Percentage of identical matches":2,"Alignment length":3,"Number of mismatches":4,"Number of gap openings":5,"Start of alignment in query":6,"End of alignment in query":7,"Start of alignment in subject":8,"End of alignment in subject":9,"Expected value":10,"Bit score":11}
    }[sys.argv[1][sys.argv[1].rindex(".")+1:]]

output = ""
for line in fraggene:
    if(line[0] == "#"): output += line + "\n"
    else: 
        data = {"Query ID":"","Subject ID":"","Percentage of identical matches":"","Alignment length":"","Number of mismatches":"","Number of gap openings":"","Start of alignment in query":"","End of alignment in query":"","Start of alignment in subject":"","End of alignment in subject":"","Expected value":"","Bit score":"","strand":""}
        line = line.split("\t")
        line[0] = line[0].split("_")
        if len(line[0]) > 4:
            for x in range(len(line[0])-4):
                line[0][x] = "_".join([line[0][x], line[0].pop(x+1)])
        data["Query ID"] = line[0][0]
        for column in filetype:
            data[column] = line[filetype[column]]
        data["Start of alignment in query"] = int(line[filetype["Start of alignment in query"]])+int(line[0][1])
        data["End of alignment in query"] = int(line[filetype["End of alignment in query"]])+int(line[0][1])
        data["strand"] = line[0][3]
        
        if(sys.argv[2][sys.argv[2].rindex(".")+1:] == "out"): output += str(data["Query ID"]) + "\t" + str(data["Subject ID"]) + "\t" + str(data["Percentage of identical matches"]) + "\t" + str(data["Alignment length"]) + "\t" + str(data["Number of mismatches"]) + "\t" + str(data["Number of gap openings"]) + "\t" + str(data["Start of alignment in query"]) + "\t" + str(data["End of alignment in query"]) + "\t" + str(data["Start of alignment in subject"]) + "\t" + str(data["End of alignment in subject"]) + "\t" + str(data["Expected value"]) + "\t" + str(data["Bit score"]) + "\n" 
        elif(sys.argv[2][sys.argv[2].rindex(".")+1:] == "bed"): output += str(data["Query ID"]) + "\t" + str(data["Start of alignment in query"]) + "\t" + str(data["End of alignment in query"]) + "\t" + str(data["Subject ID"]) + "\t" + str(data["Percentage of identical matches"]) + "\t" + str(data["strand"]) + "\n"

f = open(sys.argv[2],"w")
f.write(output[:len(output)-1])
f.close()
