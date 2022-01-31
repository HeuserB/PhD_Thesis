import sys

from PyQt6.QtWidgets import (QWidget, QCheckBox, QApplication, QVBoxLayout, QLabel)

class MultipleCheckbox(QWidget):

    def __init__(self, runs):

        super().__init__()


        # Set the label text for the user
        lb = QLabel("Select your favorite food(s):", self)
        lb.setGeometry(20, 20, 100, 20)
        lb.move(20, 20)


        # Create three checkboxes
        self.checkboxes = []
        for run_id in range(len(runs)):

            self.checkboxes.append(QCheckBox(str(runs[run_id]), self))
            cb = self.checkboxes[-1]
            cb.move(20, 70+20*run_id)
            cb.stateChanged.connect(lambda: self.Selected_Value(cb))
            
        self.label = QLabel('Nothing Selected')
        self.label.move(20, 20*len(runs))
        # Set the vertical Qt Layout

        vbox = QVBoxLayout()

        vbox.addWidget(lb)

        for run_id in range(len(runs)):
            cb = self.checkboxes[run_id]
            vbox.addWidget(cb)

        vbox.addWidget(self.label)
        self.setLayout(vbox)
        self.setWindowTitle('Form with Multiple Checkboxes')
        self.setGeometry(60, 60, 350, 200)
        self.lblText = ''


        # Display the window in the center of the screen

        #win = self.frameGeometry()
        #pos = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
        #pos = QDesktopWidget().availableGeometry().center()
        #win.moveCenter(pos)
        #self.move(win.topLeft())
        self.show()


    # Define function to read the user's input

    def Selected_Value(self, btn):
        if self.lblText != '':
            str = self.lblText
            strArray = str.split(' ,')
            self.lblText = ''
            for val in strArray:
                if btn.text() != val:
                   if self.lblText == '':
                        self.lblText = val
                   else:
                        self.lblText += ' ,' + val


            if btn.isChecked() == True:
                if self.lblText == '':
                    self.lblText = btn.text()
                else:
                    self.lblText += ' ,' + btn.text()

        else:
            if btn.isChecked() == True:
                if self.lblText == '':
                    self.lblText = btn.text()
                else:
                    self.lblText += ' ,' + btn.text()

        self.label.setText('You have selected \n' + self.lblText)
