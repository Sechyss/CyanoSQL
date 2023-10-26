import matplotlib.pyplot as plt

# Change variables and values
labels = ['Cloud', 'Shell', 'SoftCore']
values = ['1332', '1403', '59']
colours = ["#191970", "#00688B", "#528B8B"]

explode = (0.05, 0.05, 0.05)  # To highlight the Softcore fraction
fig1, ax1 = plt.subplots()
patches, texts, autotexts = ax1.pie(values, radius=2.5, explode=explode, labels=labels, rotatelabels=False,
                                    textprops={'fontsize': 14},
                                    colors=colours, autopct='%1.1f%%', shadow=False, startangle=180, counterclock=False)

for t in texts:
    t.set_color('black')
    t.set_fontsize(16)
    t.set_weight('bold')

for at in autotexts:
    at.set_color('white')
    at.set_weight('bold')

autotexts[2].set_fontsize(10)

# Donut
hole = plt.Circle((0, 0), 1, fc='white')
fig = plt.gcf()
fig.gca().add_artist(hole)

# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')
plt.tight_layout()

# Change title of the file to save the figure
plt.savefig("/Users/u1866168/Documents/OneDrive - University of Warwick/"
            "Thesis/Figures/Myoviruses_fractions_PieChart.pdf")
plt.show()
