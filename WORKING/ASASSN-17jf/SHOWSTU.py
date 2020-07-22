import matplotlib.pyplot as plt
import hipercam as hcam

data27 = hcam.hlog.Hlog.read("Quality_Reductions/2019_09_27-run015.log")
data29 = hcam.hlog.Hlog.read("Quality_Reductions/2019_09_29-run015.log")

comp_ap = '5'

gband27 = data27.tseries('2', '1') / data27.tseries('2', comp_ap)
gband29 = data29.tseries('2', '1') / data29.tseries('2', comp_ap)


fig, ax = plt.subplots()

gband27.t = gband27.t - gband27.t.mean()
gband29.t = gband29.t - gband29.t.mean()
gband27.mplot(ax, colour='red')
gband29.mplot(ax, colour='black')

plt.show()
