import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
import math
import csv
from random import random
from collections import Counter
# global circle_radius
# global obs_x_position
# global theo_x_position
# global max_x_range
# global max_y_range
# global max_lineWith
# global min_lineWith
circle_radius = 0.01
obs_x_position = 0.3
theo_x_position = 0.9
max_x_range = 1
max_y_range = 4
min_lineWith = 1
max_lineWith = 5
colorIndex = []


def drawcircl(obs, mz, w, index, ax, is_matched=0):
    if obs == 2:
        # circle1 = plt.Circle((obs_x_position, mz), circle_radius, color='r')
        circle2 = Ellipse((obs_x_position, mz), width=circle_radius,
                          height=circle_radius * 2.8, angle=0, color='r')
        ax.add_artist(circle2)
        plt.plot([obs_x_position, obs_x_position - w], [mz, mz],
                 color='k', linestyle='-', linewidth=1)
    elif obs is True:
        # circle1 = plt.Circle((theo_x_position, mz), circle_radius, color='b')
        circle2 = Ellipse((theo_x_position, mz), width=circle_radius,
                          height=circle_radius * 2.8, angle=0, color='b')
        ax.add_artist(circle2)
        if is_matched:
            ax.text(theo_x_position + 0.05, mz / 4 + 0.005, w, transform=ax.transAxes,
                    fontsize=8, verticalalignment='top', bbox=dict(boxstyle='round', facecolor=colorIndex[index], alpha=0.5), rotation = 'vertical')
        else:
            ax.text(theo_x_position + 0.05, mz / 4 + 0.005, w, transform=ax.transAxes,
                    fontsize=8, verticalalignment='top', rotation = 'vertical')

    else:
        # circle1 = plt.Circle((theo_x_position, mz), circle_radius, color='b')
        circle2 = Ellipse((theo_x_position, mz), width=circle_radius,
                          height=circle_radius * 2.8, angle=0, color='k')
        ax.add_artist(circle2)
        if is_matched:
            ax.text(theo_x_position + 0.05, mz / 4 + 0.005, w, transform=ax.transAxes,
                    fontsize=8, verticalalignment='top', bbox=dict(boxstyle='round', facecolor=colorIndex[index], alpha=0.5), rotation = 'vertical')
        else:
            ax.text(theo_x_position + 0.05, mz / 4 + 0.005, w, transform=ax.transAxes,
                    fontsize=8, verticalalignment='top', rotation = 'vertical')


def drawline(obs, theo, w, is_matched):

    if not is_matched:
        plt.plot([obs_x_position, theo_x_position], [obs, theo],
                 color='k', linestyle='dotted', linewidth=w)
    elif is_matched == 1:
        #if w == 5:
        #    plt.plot([obs_x_position, theo_x_position], [obs, theo],
        #        color='g', linestyle='-', linewidth=w)
        #else:            
        plt.plot([obs_x_position, theo_x_position], [obs, theo],
                 color='r', linestyle='-', linewidth=w)
    else:
        plt.plot([obs_x_position, theo_x_position], [obs, theo],
                 color='g', linestyle='-', linewidth=w)    


def getCordx(mz, min_cord, max_cord):
    # print mz
    # mz is list of all mzs
    if len(mz) == 0:
        return {}
    elif len(mz) == 1:
        return {mz[0]: min_cord}
    mz2 = sorted(list(set(mz)))
    r = [0] * len(mz)
    r[0] = 0
    for i in xrange(1, len(mz)):
        r[i] = mz2[i] - mz2[i - 1]
        r[i] = r[i - 1] + math.log(1 + r[i])

    curMin = r[0]
    curMax = r[-1]
    ratio = (max_cord - min_cord) / (curMax - curMin)
    return {y: ratio * (x - curMin) + min_cord for x, y in zip(r, mz2)}


def getLineWidth(w, minLineWidth, maxLineWidth):
    w2 = sorted(w)
    if len(w2) == 0:
        return {}
    elif len(w2) == 1:
        return {w2[0]: minLineWidth}
    r = [math.log(x) for x in w2]
    curMin = r[0]
    curMax = r[-1]
    ratio = (maxLineWidth - minLineWidth) / (curMax - curMin)
    return {y: ratio * (x - curMin) + minLineWidth for x, y in zip(r, w2)}


def getLineLength(w, maxLineWidth):
    w2 = sorted(w)
    if len(w2) == 0:
        return 1
    return maxLineWidth / max(w2)


def drawCurve(x0, y0, y1, is_matched):
    y = np.linspace(y0, y1, 100)
    y_mid = (y0 + y1) / 2
    y_sep = abs(y1 - y0)
    x = []
    height = 0.1 + 0.5 * y_sep / max_y_range
    # height = 0.2
    for i in y:
        t1 = (i - y_mid) * (i - y_mid) - y_sep * y_sep / 4
        t2 = y_sep * y_sep / 4
        if t2 > 0:
            t = height * t1 / t2
        else:
            return
        x.append(x0 + height * t)
    plt.plot(x, y, linewidth=0.5 + 3 * is_matched)


def drawAll(datas, num_by, filename, differ):
    # note we must use plt.subplots, not plt.subplot
    fig, ax = plt.subplots(figsize=(10, 20))

    obs_point = list(set([(x[0], x[3], x[4]) for x in datas]))
    theo_point = list(set([(x[1] / 2 * 2, x[5]) for x in datas]))
    theoid2mz = {}
    for i in theo_point:
        theoid2mz[i[0]] = i[1]
    matchedTheo2 = set(x[1] / 2 * 2 for x in datas if x[6] == 1)
    matchedTheo = set()
    for i in matchedTheo2:
        if i + 2 * num_by in matchedTheo2:
            matchedTheo.add(i)
            matchedTheo.add(i + 2 * num_by)

    print matchedTheo
    all_mz = list(set([x[1] for x in obs_point]).union(
        set([x[1] for x in theo_point])))

    mapMz2Cord = getCordx(all_mz, 2 * circle_radius,
                          max_y_range - 2 * circle_radius)

    intensityRatio = getLineLength(
        [x[2] for x in obs_point], 0.2 * max_x_range)
    for x in obs_point:
        drawcircl(2, mapMz2Cord[x[1]], x[2] * intensityRatio, 0, ax)

    # obs_pair = []
    obs_point.sort(key=lambda x: x[1])

    mass_by = 0
    for i in theo_point:
        ifbreak = 0
        for j in theo_point:
            if j[0] - i[0] == num_by * 2:
                mass_by = i[1] + j[1]
                # print i
                # print j
                ifbreak = 1
                break
        if ifbreak:
            break

    mz_begin, mz_end = 0, len(obs_point) - 1
    while mz_begin < len(obs_point) and mz_end >= 0:
        if abs(2 * obs_point[mz_begin][1] + obs_point[mz_end][1] - mass_by - 1) <= 0.1:
            temp1 = 0
            for i in matchedTheo:
                if abs(theoid2mz[i] - obs_point[mz_end][1]) < 0.1 or abs(theoid2mz[i] - obs_point[mz_begin][1]) < 0.1:
                    temp1 = 1
                    break
            drawCurve(0.3, mapMz2Cord[obs_point[mz_begin][1]], mapMz2Cord[obs_point[mz_end][1]], temp1)
            mz_begin += 1
            mz_end -= 1
        elif 2 * obs_point[mz_begin][1] + obs_point[mz_end][1] - mass_by - 1 <= 0:
            mz_begin += 1
        else:
            mz_end -= 1

    for x in theo_point:
        # for x in zip(theo_point, label):

        if x[0] < num_by * 2:
            drawcircl(x[0] < 2 * num_by, mapMz2Cord[x[1]],
                      "b%d" % (x[0] / 2 + 1), x[0] / 2, ax, x[0] in matchedTheo)
        else:
            drawcircl(x[0] < 2 * num_by, mapMz2Cord[x[1]],
                      "y%d" % (num_by - x[0] / 2 + num_by), x[0] / 2 - num_by, ax, x[0] in matchedTheo)
    mapW2lineWide = getLineWidth(
        list(set([x[2] for x in datas])), min_lineWith, max_lineWith)

    for x in datas:
        drawline(mapMz2Cord[x[3]], mapMz2Cord[x[5]], mapW2lineWide[x[2]], x[6])
    #for i in differ:
    #    print i
    #    drawline(mapMz2Cord[datas[i][3]], mapMz2Cord[datas[i][5]], mapW2lineWide[datas[i][2]], 2)    

    ax.set_xlim([0, max_x_range])
    ax.set_ylim([0, max_y_range])
    ax.text(obs_x_position - 0.01, 1.1, "Observed, V", transform=ax.transAxes,
            fontsize=14, verticalalignment='top', rotation = 'vertical')
    ax.text(theo_x_position - 0.01, 1.1, "Theoretical, U", transform=ax.transAxes,
            fontsize=14, verticalalignment='top', rotation = 'vertical' )

    # drawcircl(0, 10, ax)

    plt.axis('off')
    # drawCurve(0.5, 0, 0.1)
    plt.savefig(filename)
    plt.clf()
    plt.close()
    # fig = matplotlib.pyplot.gcf()
    # fig.set_size_inches(18.5, 10.5)
    # ig.savefig('test2png.png', dpi=100)


def read_and_process_one_dataset(inFile, outFile):
    dataset = []
    num_bion = []
    sids = []
    with open(inFile, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            if int(row[0]) < 0:
                dataset.append([])
                num_bion.append(-int(row[0]))
                sids.append(int(row[1]))
                continue
            dataset[-1].append((int(row[0]), int(row[1]), float(row[2]),
                                float(row[3]), float(row[4]), float(row[5]), int(row[6])))
    # match = set([4, 9])
    # unmatch = set()

    # for i in xrange(len(dataset[19])):
    #     if dataset[19][i][1] / 2 in match or dataset[19][i][1] / 2 - num_bion[19] / 4 in match:
    #         dataset[19][i] = (dataset[19][i][0], dataset[19][i][1], dataset[19][i][
    #                           2], dataset[19][i][3], dataset[19][i][4], dataset[19][i][5], 1)
    #     if dataset[19][i][1] / 2 in unmatch or dataset[19][i][1] / 2 - num_bion[19] / 4 in unmatch:
    #        dataset[19][i] = (dataset[19][i][0], dataset[19][i][1], dataset[19][i][
    #                           2], dataset[19][i][3], dataset[19][i][4], dataset[19][i][5], 0)
    differ1 = []
    differ2 = []
    for i in xrange(len(dataset[0])):
        if not dataset[0][i][-1] == dataset[1][i][-1]:
            print i
            if dataset[0][i][-1] == 1:
                differ1.append(i)
            else:
                differ2.append(i)    
    print differ1, differ2            
    
    #drawAll(dataset[0], num_bion[0] / 4, '%s/%d.pdf' % (outFile, 0), differ1)
    #drawAll(dataset[1], num_bion[1] / 4, '%s/%d.pdf' % (outFile, 1), differ2)
    #return dataset
    maxDraw = 10
    if maxDraw >= 0:
        maxDraw = min(maxDraw, len(dataset))
        for i in xrange(maxDraw):
        	drawAll(dataset[i], num_bion[i] / 4, '%s/%d.pdf' % (outFile, i), [])
            #if not i == 19:
            #    continue
            #if sids[i] in diff_sid:
            #    print "drawing %d" % i
            
    return dataset




for i in xrange(100):
    colorIndex.append((random(), random(), random()))

dataset1 = read_and_process_one_dataset('matchings_one.txt', 'figures')

# dataset2 = read_and_process_one_dataset('matchings1.txt', 'figures0')
# dataset2 = read_and_process_one_dataset('matchings0.txt', 'figure')

diff = 0



# 57577.6 59123.4
# 47198.2 48781.3
