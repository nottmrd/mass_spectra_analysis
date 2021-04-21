# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:54:32 2020

@author: peter
"""

"""
This is a tool developed for my father, a Chemistry PhD, to guess at the chemical formula of a 
certain molecule using Mass Spectra data (M+ 1 peaks).

He knew the weight of the molecule, and that it had a Benzyne ring.
Also he knew that there would be more Hydrogen atoms than Carbon, and that the molecule
would be composed of (at most) Carbon, Hydrogen, Oxygen, and Chlorine.
And further, that the element being in the Alkyl group would have a formula of C(n)H(2n+1).

In addition to this, he knew that a molecule broke off at a certain point, and the molecule that
broke off would have to follow the same constraints as above, it just would be smaller (i.e. less
Carbon and Hydrgoen atoms).

With these constraints, this program formulates what possibilites the molecule could be.
"""


import itertools

#After all possibilites are calculated, this function sorts the list by Variance (the difference 
#between the target weight and what is calculated). We consider these to be the most likely answers.
def mySort(list):
    return list["Variance"]

#This function takes into account all 4 elements. The other assumes it's likely there isn't Chlorine.
def weight_cl(target, left, right, righter, next_peak, best_answer=0):
    ratio = right/left
    ratio_cl = righter/left
    #A guess at the maximum amount of Hydrogen atoms, very unlikely to be more than this many
    max_h = 20
    #The proportion of these elements that have an extra electron (M+1, M+2)
    c_ratio = 1.11
    h_ratio = 0.015
    o_ratio = 0.04
    cl_ratio = 25
    #How many electrons each element needs to make the happy stable molecules. See the periodic table 
    c_spare = 4
    h_spare = 0
    o_spare = 2
    cl_spare = 1
    #Atomic weight
    chemicals = [12, 1, 16, 35]
    #We find all the combinations using the atomic weights, that could sum up to our target weight
    combos = [seq for i in range(100, 0, -1) for seq in itertools.combinations_with_replacement(chemicals, i) if sum(seq) == target]
    answers = []
    #We perform further analysis on these comboniations and determine their makeup.
    for each in combos:
        c = 0
        h = 0
        o = 0
        cl = 0
        for i in each:
            if i == 1:
                h += 1
            elif i == 12:
                c += 1
            elif i == 16:
                o += 1
            else:
                cl += 1
        if h <= max_h:
            benzyne = False
            c_ans = (c_ratio/(100+c_ratio)) * c
            h_ans = (h_ratio/(100+h_ratio)) * h
            o_ans = (o_ratio/101.24) * o
            cl_ans = (cl_ratio/100) * cl
            if c >= 6 and h >= 4:
                benzyne = True
            ans = c_ans + h_ans + o_ans
            var_cl = round(ratio_cl - cl_ans, 3)
            var = round(ratio - ans,3)
            if ((c_spare * c)+(h_spare * h)+(o_spare * o)+(cl_spare * cl)) % 8 == 0:
                balanced = True
            else:
                balanced = False
            if h == (c*2)+1:
                cnh2n1 = True
            else:
                cnh2n1 = False
            if h > c:
                more = True
            else:
                more = False
            if benzyne == True and balanced == True and more == True and o <= 4:
                answers.append({"C":c, "H":h, "O":o, "Cl":cl, "C(n)H(2n+1)":cnh2n1, "Variance":var, "Variance M+2":var_cl, "Benzyne ring":benzyne, "Balanced":balanced, "More H than C":more})
    answers.sort(key=mySort)
    print("")
    print("Combination of chemical elements for a total weight of:", target)
    print("Constrained by a maximum number of hydrogens of:", max_h)
    print("")
    #There is a second peak that indicates some of the molecule broke off. 
    #This function formulates quite simply what this breakoff molecule could be. 
    #A glance at the answers this function produces will tell if these break off molecules are even possible.
    for each in answers:
        #Memo is used to account for duplicates. We don't need duplicate chemical formulas.
        memo = []
        new = target-next_peak
        combos = [seq for i in range(100, 0, -1) for seq in itertools.combinations_with_replacement(chemicals, i) if sum(seq) == new]
        for each in combos:
            c = 0
            h = 0
            o = 0
            cl = 0
            for i in each:
                if i == 12:
                    c += 1
                elif i == 1:
                    h += 1
                elif i == 16:
                    o += 1
                elif cl == 35:
                    cl += 1
            #We make sure that this possible molecule that broke off is possible (i.e. there are enough Carbon atoms
            #in one of our best case answers, that a breakoff of these carbon atoms would be possible.)
            if c <= answers[best_answer]["C"] and h <= answers[best_answer]["H"] and o <= answers[best_answer]["O"] and cl <= answers[best_answer]["Cl"]:
                if ("C:",c,"H:",h,"O:",o,"Cl:",cl) not in memo:
                    memo.append(("C:",c,"H:",h,"O:",o,"Cl:",cl))
    for each in answers:
        print("   ", each)
    print("")
    print("Breakdown before the peak at:",next_peak,". Returns these options:")
    print("")
    print("   ", memo)
    
#This function assumes it's unlikely that there is any Chlorine.
def weight(target, left, right, next_peak, best_answer=0):
    ratio = right/left
    #A guess at the maximum amount of Hydrogen atoms, very unlikely to be more than this many
    max_h = 50
    #The proportion of these elements that have an extra electron (M+1, M+2)
    c_ratio = 1.11
    h_ratio = 0.015
    o_ratio = 0.04
    #How many electrons each element needs to make the happy stable molecules. See the periodic table 
    c_spare = 4
    h_spare = 0
    o_spare = 2
    #Atomic weight
    chemicals = [12, 1, 16]
    #We find all the combinations using the atomic weights, that could sum up to our target weight
    combos = [seq for i in range(100, 0, -1) for seq in itertools.combinations_with_replacement(chemicals, i) if sum(seq) == target]
    answers = []
    #We perform further analysis on these comboniations and determine their makeup.
    for each in combos:
        c = 0
        h = 0
        o = 0
        for i in each:
            if i == 1:
                h += 1
            elif i == 12:
                c += 1
            else:
                o += 1
        if h <= max_h:
            benzyne = False
            c_ans = (c_ratio/(100+c_ratio)) * c
            h_ans = (h_ratio/(100+h_ratio)) * h
            o_ans = (o_ratio/(101.24)) * o
            if c >= 6 and h >= 4:
                benzyne = True
            ans = c_ans + h_ans + o_ans
            var = round(ratio - ans,3)
            if ((c_spare * c)+(h_spare * h)+(o_spare * o)) % 8 == 0:
                balanced = True
            else:
                balanced = False
            if h == (c*2)+1:
                cnh2n1 = True
            else:
                cnh2n1 = False
            if h > c:
                more = True
            else:
                more = False
            if benzyne == True and balanced == True and more == True and o <= 4:
                answers.append({"C":c, "H":h, "O":o, "C(n)H(2n+1)":cnh2n1, "Variance":var, "Benzyne ring":benzyne, "Balanced":balanced, "More H than C":more})
    answers.sort(key=mySort)
    print("")
    print("Combination of chemical elements for a total weight of:", target)
    print("Constrained by a maximum number of hydrogens of:", max_h)
    print("")
    #There is a second peak that indicates some of the molecule broke off. 
    #This function formulates quite simply what this breakoff molecule could be. 
    #A glance at the answers this function produces will tell if these break off molecules are even possible.
    for each in answers:
        #Memo is used to account for duplicates. We don't need duplicate chemical formulas.
        memo = []
        new = target-next_peak
        combos = [seq for i in range(100, 0, -1) for seq in itertools.combinations_with_replacement(chemicals, i) if sum(seq) == new]
        for each in combos:
            c = 0
            h = 0
            o = 0
            for i in each:
                if i == 12:
                    c += 1
                elif i == 1:
                    h += 1
                elif i == 16:
                    o += 1
            #We make sure that this possible molecule that broke off is possible (i.e. there are enough Carbon atoms
            #in one of our best case answers, that a breakoff of these carbon atoms would be possible.)
            if c <= answers[best_answer]["C"] and h <= answers[best_answer]["H"] and o <= answers[best_answer]["O"]:
                if ("C:",c,"H:",h,"O:",o) not in memo:
                    memo.append(("C:",c,"H:",h,"O:",o))
    for each in answers:
        print("   ", each)
    print("")
    print("Breakdown before the peak at:",next_peak,". Returns these options:")
    print("")
    print("   ", memo)


#These are actual Mass Spectra results.

# weight(176, 218, 31, 161, 0)
weight(206, 21, 3, 161, 1)
# weight(320, 111, 27)
# weight(320, 4, 1)
# weight_cl(196, 96, 14, 35, 161, 0)
    