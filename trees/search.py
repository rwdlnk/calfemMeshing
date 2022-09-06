
pts = [[1,6,7,2],
       [2,7,8,3],
       [3,8,9,4],
       [4,9,10,5],
       [9,12,11,10],
       [8,13,12,9],
       [7,14,13,8],
       [6,15,14,7],
       [15,16,17,14],
       [14,17,18,13],
       [13,18,19,12],
       [12,19,20,11]]
item = 15
try:
  #search for the item
  index = pts[5].index(item)
  print('The index of', item, 'in the list is:', index)
except ValueError:
  print('item not present')
