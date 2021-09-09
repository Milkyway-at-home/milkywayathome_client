import tkinter as tk
import os
import re
import matplotlib
from matplotlib.figure import Figure
import numpy as np
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
#-------------------------------------------------------------------------------
#Created by: Tom Donlon, Shreyas Seethalla, Rene Saenz
#read in settings from lua file
#and construct a dictionary to access them
#-------------------------------------------------------------------------------
#bad way to do the filename stuff, but I couldn't get os to play nice
dirname = os.path.dirname(os.path.abspath(__file__))
dirname = dirname[:dirname.rindex('/')] #each time removes last directory from dirname
dirname = dirname[:dirname.rindex('/')]
filename = dirname + '/nbody/sample_workunits/settings2.lua'
filename2= filename #grabs filename for input
filename2 = filename2[:filename2.rindex('/')]
with open(filename, 'r') as file:
    data = file.readlines()

#build dict with index of line and ach relevant setting
settings = {}

for i in range(len(data)):
    line = data[i]
    
    if line != '\n':
        if (line.split() and line.split()[0] == 'function'): #first line that isn't commented starting
                                          #with actual lua code is where settings stop
            break
        elif line.lstrip()[:2] != '--':
            settings[i] = line
            
#Splicing form:
##nbodyMinVersion       |~= ~|"1.80"|  | ~--~| |~--~| MINIMUM APP VERSION |~$~ |entry | 1.80 |~^~ |4 |~*~ |.5
#Tildas surround what is not included in dictionary
#Splits before equal sign, around value, around the comment, interface type, and limits

#the keys in the dict are the original line numbers that each setting comes from
for key in settings:
    tmp_list1 = settings[key].strip().split('--')
    tmp_list2 = tmp_list1[0].split('=')
    tempStr = tmp_list2[-1]
    tempStr = tempStr[1:]
    tmp_list2[-1] = tempStr
    tmp_list3 = tmp_list2[-1].split(' ', 1) # Used to isolate the value
    tmp_list4 = re.split(r"[$*^|]\s*",tmp_list1[-1])
    settings[key] = tmp_list2[0:1] + tmp_list3[:] + tmp_list1[1:2] + tmp_list4[:]

#0 - name w/ spaces
#1 - value
#2 - spaces
#3 - spaces
#4 - comment
#5 - interface type
#6 - default value
#7 - maximum value
#8 - minimum value





#BreakPoint makes a (sub)title. Called breakpoint because both
#in settings.lua and in this GUI they are breaks from settings.
def BreakPoint(name, position, frameName):
    aLabel = tk.Label(frameName, text=name)
    aLabel.pack(row=position, column = 0)
    
#Textual entry boxes. All update when return is pressed.
def EntryBuild(name, position, current_default, upper, downer, current_key, box_width, frameName):
    current_default = current_default.strip()
    def update(self):
        nonlocal a
        settings[current_key][1] = a.get()
        a.config(text = settings[current_key][1])
    defaultString = tk.StringVar(root, value = current_default)
    a = tk.Entry(frameName, width=box_width, textvariable= defaultString)
    aLabel = tk.Label(frameName, text= name)
    a.grid(row=position, column = 2)
    aLabel.grid(row = position, column = 0)
    #Creates a limit label to let users know what can be inputed
    limits = f"Min: {downer}; Max: {upper.strip()}"
    aLimit = tk.Label(frameName, text = limits)
    aLimit.grid(row = position, column = 4, padx=10) #Adds limits to the left of the entries
    root.bind("<Return>", update)

    
#Clickable buttons, update the dataframe on toggle
def ButtonBuild(name,  position,  current_default, current_key, frameName):
    current_default = 'true' if current_default == '1' else 'false' #lmits use binary, need to convert
    def aToggle():
        nonlocal aButton
        # remember to figure out the correct syntax for this (and comment your code correctly)
        if aButton.config('text')[-1] == 'true':
            aButton.config(text='false')
            settings[key][1] = "false"
        else:
            aButton.config(text='true')
            settings[key][1] = "true"
    
    boolean = current_default
    global settings
    
    aButton = tk.Button(frameName, text = current_default, command = aToggle)
    aLabel = tk.Label(frameName, text=name)
    aButton.grid(row = position, column = 2)
    aLabel.grid(row = position, column = 0)
    
#Increment/Decrement readout
def ValuesBuild(name, position, default,  upper, downer, current_key, frameName):
    
    def Increment(position, value, higher, lower, current_key): # also remember to use that same correct syntax here.
        
        global settings
        nonlocal aValue
        nonlocal tempval
        
        if (tempval< higher):
            tempval+= 1
        aValue.config(text=tempval)
        
        
        default = tempval
        settings[current_key][1]= tempval
        return tempval
        
    def Decrement(position, value, higher, lower, current_key):
        
        global settings
        nonlocal aValue
        nonlocal tempval

        if (tempval > lower):
            tempval -= 1
        aValue.config(text= tempval)
        default = tempval
        settings[current_key][1] = tempval
        
        return tempval
    
    tempval = int(default)
    aTitle = tk.Label(frameName, text = name)
    aValue = tk.Label(frameName, text = default)
    aIncrement =tk.Button(frameName, text = ">", command = lambda:Increment(position, int(tempval), int(upper), int(downer), key))
    aDecrement= tk.Button(frameName, text = "<", command = lambda:Decrement(position, int(tempval), int(upper), int(downer), key))
    aValue.pack(row = position, column = 2)
    aTitle.pack(row=position, column = 0)
    aIncrement.pack(row = position, column = 3)
    aDecrement.pack(row = position, column = 1)
    #Creates a limit label to let users know what can be inputed
    limits = f"Min: {downer}; Max: {upper.strip()}"
    aLimit = tk.Label(text = limits)
    aLimit.pack(row = position, column = 4, padx=10) #Adds limits to the left of the entries
    

def onFrameConfigure(canvas):
    '''Reset the scroll region to encompass the inner frame'''
    canvas.configure( scrollregion=canvas.bbox("all"))

root = tk.Tk()
canvas = tk.Canvas(root, borderwidth=0, background="light gray")
frame = tk.Frame(canvas, background="light gray")
vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview) #creates scrollbar to far right
canvas.configure(yscrollcommand=vsb.set)

vsb.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)
canvas.create_window((4,4), window=frame, anchor="nw")

frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

pos = 0 #represents which row widgets appear on

#Creates the widgets in the frame based on the interface type specified
for key in settings:
    
    if (len(settings[key])<9):
        continue
    
    if settings[key][5] == "break ":
        BreakPoint(settings[key][0], pos, frame)
        pos = int((pos/100 +1)*100)
    if settings[key][5] == "button ":
        ButtonBuild(settings[key][0], pos, settings[key][6].strip(), key, frame)
        pos += 1
    if settings[key][5] == "short ":
        ValuesBuild(settings[key][0], pos, settings[key][6], settings[key][7], settings[key][8], key, frame)
        pos += 1
        continue
    if settings[key][5] == "l-entry " or settings[key][5] == "l-q-entry ":
        EntryBuild(settings[key][0], pos, settings[key][6], settings[key][7],  settings[key][8], key, 12, frame)
        pos+=1
    if settings[key][5] == "entry " or settings[key][5] == "q-entry ":
        EntryBuild(settings[key][0], pos, settings[key][6], settings[key][7],  settings[key][8], key, 6, frame)
        pos+=1
        
childList = canvas.winfo_children #creates a list of all widgets in the first (most shallow) level

fig = Figure(figsize=(5, 4), dpi=100)
t = np.arange(0, 3, .01)
ax =fig.add_subplot(projection='3d')
filename2 += '/manual_bodies_example.in'
with open(filename2, 'r') as file:
  lines=file.readlines()
coords = []

for i in range(len(lines)):
  if(i>0):
    coords.append([float(x) for x in lines[i].split()[2:]])

x, y, z = [i[0] for i in coords],[i[1] for i in coords],[i[2] for i in coords]


u,v,w = [i[3] for i in coords],[i[4] for i in coords],[i[5] for i in coords]
ax.quiver(x, y, z, u, v, w, length = 20,normalize = "True")


canvas2 = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas2.draw()
canvas2.get_tk_widget().pack()
root.wm_title("NBody Lite")

#-------------------------------------------------------------------------------
#print out settings to file when you click the run button
#-------------------------------------------------------------------------------

def updateSettings():
    def readFromSettings(i):
        if(settings[i][5].strip() == 'q-entry' or settings[i][5].strip() == 'l-q-entry'): #Adds quotes around designated settings
            settings[i][1] = '\"' + settings[i][1].strip() + '\"'

        file.write(f'{settings[i][0]}= {settings[i][1]}{settings[i][2]} -- --{settings[i][4]}$ {settings[i][5]}| {settings[i][6]}^ {settings[i][7]}* {settings[i][8]}\n')

    filename2 = dirname + '/nbody/sample_workunits/settings2.lua'
    print("Updating settings...")
    with open(filename2, 'w') as file:
        endIndex = 0
        for i in range(0, len(data)):
            line = data[i]
            if line == '\n':
                file.write('\n')
            elif line[0] == '-' or line[0] == ' ' or not line: #line is a comment
                file.write(line)
            elif line.split()[0] == 'function': #Ends the loop and copies rest from old lua file
                endIndex = i
                break
            else:
                readFromSettings(i)


        for i in range(endIndex, len(data)):
            file.write(data[i])


def retrieve_input(childList):

    #Remove label objects
    childList = [child for child in childList if child.winfo_class() != 'Label']

    #iterate through entries
    for key, childIndex in zip(settings, range(0, len(childList))):
        child = childList[childIndex]
        if child.cget('text') == '<': #for short interface, not in use
            childIndex -= 1
            continue
        if child.cget('text') == '>': #for short interface, not in use
            row = child.pack_info()['row']
            settings[key][1] = window.pack_slaves(row, 2)[0].cget('text')
            continue
        if child.winfo_class() == 'Entry':
            settings[key][1] = child.get()
        if child.winfo_class() == 'Button':
            settings[key][1] = child.cget('text')


    updateSettings()

buttonCommit=tk.Button(root, height=1, width=10, text="Apply",
                    command=lambda: retrieve_input(childList[0].winfo_children()))
buttonCommit.pack(side="bottom")
root.mainloop()
