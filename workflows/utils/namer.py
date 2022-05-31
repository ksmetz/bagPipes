def namer(sampleTable, columnNames):
    nameList = []
    # For columns which you want to list ONLY if there is a single one
    for column in columnNames:
        columnSet = list(set(sampleTable[column]))
        if len(columnSet) == 1:
            nameList.append(columnSet[0])
    return("_".join(nameList))