"""Script to test APO Exposure timer calc"""

from APOExptime import Sky, Target, Instrument, Observation, makeplots

sky = Sky(lunar_phase=0, seeing=1)
star = Target(21, 'VEGAMAG', [5000, 6000], temp=30000)
inst1 = Instrument('Arctic')


ob1 = Observation(star, sky, inst1)


sn1 = ob1.SNfromTime(100)
print(sn1)

t1 = ob1.TimefromSN(100)

#plt.plot(t1[0], t1[1])
# plt.plot(times[0], times[1])

#plt.plot(sn1[2][0], sn1[2][1])
# plt.plot(times[0], times[1])


makeplots(ob1, 'SN')
