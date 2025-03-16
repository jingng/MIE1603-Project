import matplotlib.pyplot as plt


def read_file(filepath):
  with open(filepath, 'r') as file:
    data = []
    for line in file:
      newline=line.split()
      data.append(newline)

  return data

def read_values(data):

  sol_dict={'t':[], 'bound':[],'obj':[], 'gap':[]}

  for i in range(len(data)):

    time = float(data[i][0].split(',')[0])
    sol_dict['t'].append(time)

    obj = float(data[i][2])
    sol_dict['obj'].append(obj)

    bound = float(data[i][1].split(',')[0])

    if bound > 0.05:
        sol_dict['bound'].append(bound)
        gap = abs(bound - obj)/obj
        sol_dict['gap'].append(gap)

  return sol_dict

# filepath = '/Users/jingong/Github/MIE1603-Project/experiment_jess_3600/'
filepath = 'C://Users//Jin//Github//MIE1603-Project//experiment_jess_3600//'

instances = [
# #   'instance-jess-min_1000_test.csv',
#   'instance-jess-min_5000_test.csv',
#   'instance-jess-min_8000_test.csv',
#   'instance-jess-min_10000_test.csv',
#   'instance-jess-min_15000_test.csv'

  #   'instance-jess-min_1000_test.csv',
  'instance-jess-min-rem-aggressive_5000.csv',
  'instance-jess-min-rem-aggressive_8000.csv',
  'instance-jess-min-rem-aggressive_10000.csv',
  'instance-jess-min-rem-aggressive_15000.csv'
]
lengths = ['5 km','8 km','10 km','15 km']

for i in range(len(instances)):

    data = read_file(filepath+instances[i])
    data_dict = read_values(data)

    time_first_bound = len(data_dict['t']) - len(data_dict['gap'])
    x = data_dict['t'][time_first_bound:]
    bound = data_dict['bound']
    obj = data_dict['obj'][time_first_bound:]
    gap = data_dict['gap']

    plt.plot(x, gap, label=lengths[i])
    plt.legend()
    plt.title('MIP Gap for North Delta')
    plt.xlabel('Time (s)')
    plt.ylabel('MIP Gap')
    plt.ylim([0,0.4])

plt.show()

