import itertools
from collections import Counter
from math import floor, ceil, cos, sin, pi, sqrt

import numpy


class Encoding:
    def GlobalEncoding(self, seq):

        Feat = []

        H = numpy.zeros((20, len(seq))).astype(int)


        L = 6  # size of the patch used for extracting the descriptors
        if len(seq) < L:
            seq[len(seq): L] = seq[len(seq)]  # # to avoid proteins too short

        Len = len(seq)
        for i in range(Len):
            t = 0
            if seq[i] == 'A':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'R':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'N':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'D':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'C':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'Q':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'E':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'G':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'H':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'I':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'L':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'K':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'M':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'F':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'P':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'S':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'T':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'W':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'Y':
                H[t, i] = 1
            else:
                H[t, i] = 0

            t = t + 1
            if seq[i] == 'V':
                H[t, i] = 1
            else:
                H[t, i] = 0

        F = numpy.zeros((20, L * 2))  # for 0 and 1 for each subsequences

        # 揷omposition, " i.e., the frequency of 0s and 1s
        for i in range(20):  # the 10 binarysequence H

            S = max([len(seq) / L, 1])
            t = 0
            for j in range(1, L + 1):

                F[i, t] = round(list(H[i, :])[floor((j - 1) * S): floor(j * S)].count(1) / S, 4)
                t = t + 1
                # F[i, t] = round(list(H[i,:])[ floor((j - 1)  * S): floor((j) * S)-1].count(0)/S,4)

                if j == 1:
                    F[i, t] = round(list(H[i, :])[floor((j - 1) * S): floor((j) * S)].count(0) / S, 4)
                else:
                    F[i, t] = round(list(H[i, :])[floor((j - 1) * S): floor((j) * S) - 1].count(0) / S, 4)
                t = t + 1

        # transition? i.e., the percent of frequency with which 1 is followed by 0 or 0 is followed by 1 in a characteristic sequence

        F1 = [0] * 20  # for 0-1 transition, 1 , 11, 111 fopr each subsequences

        for i in range(20):  # the 10 binarysequence H
            S = max([len(seq) / L, 1])
            t = 0

            temp = []
            for j in range(1, L + 1):
                Sezione = list(H[i, :])[floor((j - 1) * S): floor((j * S) - 1)]
                Sezione1 = Sezione[1:len(Sezione)]
                Sezione2 = Sezione[2:len(Sezione)]

                counter1 = 0
                counter2 = 0
                counter3 = 0
                for k in range(len(Sezione1)):
                    if Sezione[k] == 1 and Sezione1[k] == 0:
                        counter1 += 1
                    if Sezione[k] == 0 and Sezione1[k] == 1:
                        counter1 += 1

                for k in range(len(Sezione1)):

                    if Sezione[k] == 1 and Sezione1[k] == 1:
                        counter2 += 1

                for k in range(len(Sezione2)):

                    if Sezione[k] == 1 and Sezione1[k] == 1 and Sezione2[k] == 1:
                        counter3 += 1
                temp.extend([counter1, counter2, counter3])

            F1[i] = temp

        F1 = numpy.array(F1).astype("float")

        F = list(F.astype("float").flatten("F"))
        F1 = list(F1.flatten("F"))
        F.extend(F1)

        return numpy.array(F)
    def Sline(self, seq):

        V = []
        C1 = 'AVLIMC'
        C2 = 'FWYH'
        C3 = 'STNQ'
        C4 = 'KR'
        C5 = 'DE'
        C6 = 'GP'

        ############################################group1-------------------

        length = len(seq)
        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C2 or seq[i] in C3:
                seq1[i] = ord('B')
            elif seq[i] in C4:
                seq1[i] = ord('J')
            elif seq[i] in C5:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[0]) ** 2 for s in seqx]
        b1 = [(s - V[1]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group2-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C2 or seq[i] in C4:
                seq1[i] = ord('B')
            elif seq[i] in C3:
                seq1[i] = ord('J')
            elif seq[i] in C5:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[4]) ** 2 for s in seqx]
        b1 = [(s - V[5]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group3-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C2 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C3:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[8]) ** 2 for s in seqx]
        b1 = [(s - V[9]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group4-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C2 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C3:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[12]) ** 2 for s in seqx]
        b1 = [(s - V[13]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group5-------------------

        length = len(seq)
        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C3 or seq[i] in C4:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C5:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[16]) ** 2 for s in seqx]
        b1 = [(s - V[17]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group6-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C3 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[20]) ** 2 for s in seqx]
        b1 = [(s - V[21]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################ group7-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C3 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[24]) ** 2 for s in seqx]
        b1 = [(s - V[25]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group8-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C4 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[28]) ** 2 for s in seqx]
        b1 = [(s - V[29]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group9-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C4 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')


        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[32]) ** 2 for s in seqx]
        b1 = [(s - V[33]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group10-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C1 or seq[i] in C5 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C2:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C4:
                seq1[i] = ord('U')

        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)

        seqx = [0] * length

        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)

        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        # ############################################group11-------------------

        length = len(seq)
        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C3 or seq[i] in C4:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C5:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group12-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C3 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group13-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C3 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C4:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group14-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C4 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group15-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C4 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group16-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C2 or seq[i] in C5 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C3:
                seq1[i] = ord('O')
            elif seq[i] in C4:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group17-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C3 or seq[i] in C4 or seq[i] in C5:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C2:
                seq1[i] = ord('O')
            elif seq[i] in C6:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group18-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C3 or seq[i] in C4 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C2:
                seq1[i] = ord('O')
            elif seq[i] in C5:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group19-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C3 or seq[i] in C5 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C2:
                seq1[i] = ord('O')
            elif seq[i] in C4:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        ############################################group20-------------------

        length = len(seq)

        seq1 = numpy.array([0] * length)
        for i in range(length):
            if seq[i] in C4 or seq[i] in C5 or seq[i] in C6:
                seq1[i] = ord('B')
            elif seq[i] in C1:
                seq1[i] = ord('J')
            elif seq[i] in C2:
                seq1[i] = ord('O')
            elif seq[i] in C3:
                seq1[i] = ord('U')
        idx1 = numpy.where(seq1 == ord("B"))[0]
        idx2 = numpy.where(seq1 == ord("J"))[0]
        idx3 = numpy.where(seq1 == ord("O"))[0]
        idx4 = numpy.where(seq1 == ord("U"))[0]
        AN = len(idx1)
        BN = len(idx2)
        CN = len(idx3)
        DN = len(idx4)
        seqx = [0] * length
        seqy = [0] * length
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0

        for i in range(length):
            if seq1[i] == ord('B'):
                n1 = n1 + 1
                seqx[i] = round(float(cos((pi / 2) * (n1 / (AN + 1)))), 4)
                seqy[i] = round(float(sin((pi / 2) * (n1 / (AN + 1)))), 4)
            elif seq1[i] == ord('J'):
                n2 = n2 + 1
                seqx[i] = round(float(cos(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
                seqy[i] = round(float(sin(pi / 2 + (pi / 2) * (n2 / (BN + 1)))), 4)
            elif seq1[i] == ord('O'):
                n3 = n3 + 1
                seqx[i] = round(float(cos(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
                seqy[i] = round(float(sin(pi + (pi / 2) * (n3 / (CN + 1)))), 4)
            elif seq1[i] == ord('U'):
                n4 = n4 + 1
                seqx[i] = round(float(cos(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
                seqy[i] = round(float(sin(3 * pi / 2 + (pi / 2) * (n4 / (DN + 1)))), 4)
        V.append(round(1 / length * numpy.sum(seqx), 4))
        V.append(round(1 / length * numpy.sum(seqy), 4))
        a1 = [(s - V[36]) ** 2 for s in seqx]
        b1 = [(s - V[37]) ** 2 for s in seqy]
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(a1)), 4))
        V.append(round(sqrt(1 / (length - 1) * numpy.sum(b1)), 4))

        return V

