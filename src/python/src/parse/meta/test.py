from knownKeywords import knownKeywords as kw

for i in kw:
  print i
  for j in kw[i]:
    print '  '+j+': '+str(kw[i][j].req)
