#!/usr/bin/env python

import pygtk
pygtk.require('2.0')
import gtk

import formatConvertor
import formatType
import gui.treefilebrowser


class PicosGtk:

    def __init__(self):
        self.init_window()
        self.create_widgets()
        self.pack_widgets()
        self.show_all()
    
    ###################   Init Subfunctions  ####################################  
    def init_window(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_title("Picos")
        self.window.connect("delete_event", self.delete_event)
        self.window.set_border_width(10)
        self.window.resize(800, 600)

        
    def create_widgets(self):
        self.button_convert = gtk.Button("Convert")
        self.button_convert.connect("clicked", self.convert, "convert")
        
        self.text_input = gtk.Entry()
        self.text_output = gtk.Entry()
        
        self.button_input = gtk.Button("Open file...")
        self.button_input.connect("clicked", self.openDialog, self.text_input)
        self.button_output = gtk.Button("Open file...")
        self.button_output.connect("clicked", self.openDialog, self.text_output)
        
        self.check_input_open = gtk.CheckButton("Half-open")
        self.check_output_open = gtk.CheckButton("Half-open")
        
        self.radio_dir = gtk.RadioButton(group=None, label="Convert directory")
        self.radio_dir.set_active(True)
        self.radio_file = gtk.RadioButton(group=self.radio_dir, label="Convert file")
        
        list_input = [formatType.ELAND, formatType.RAWELAND, formatType.BED]
        self.combo_input = gtk.combo_box_new_text()
        for item in list_input:
            self.combo_input.append_text(item)
        self.combo_input.set_active(0)
        
        list_output = [formatType.BED, formatType.WIG, formatType.MONOWIG, formatType.PEAK, formatType.SPEAK]
        self.combo_output = gtk.combo_box_new_text()
        for item in list_output:
            self.combo_output.append_text(item)
        self.combo_output.set_active(0)  
        
        # Init TreeFileBrowser
        self.file_browser = gui.treefilebrowser.TreeFileBrowser()
        self.file_browser_view = self.file_browser.get_view()
        self.file_browser_scrolled = self.file_browser.get_scrolled()
        signal = self.file_browser.connect("cursor-changed", self.dir_selected)
        
    def pack_widgets(self):
        self.Hbox1 = gtk.HBox(False, 0)
        self.Vbox1 = gtk.VBox(False, 0)
        self.box1 = gtk.HBox(False, 0)
        self.box2 = gtk.HBox(False, 0)
        self.box3 = gtk.HBox(False, 0)
        self.box4 = gtk.VBox(False, 0)
        self.window.add(self.Hbox1)
        self.separator = gtk.HSeparator()
        self.separator2 = gtk.HSeparator()
        
        self.Vbox1.add(self.box1)
        self.Vbox1.add(self.separator)
        self.Vbox1.add(self.box2)
        self.Vbox1.add(self.separator2)
        self.Vbox1.add(self.box3)
        self.Vbox1.add(self.box4)
        
        self.Hbox1.pack_start(self.file_browser_scrolled, expand=True, fill=True, padding=0)
        self.Hbox1.add(self.Vbox1)
    
        self.box1.pack_start(self.combo_input, expand=True, fill=False, padding=0)
        self.box1.pack_start(self.combo_output, expand=True, fill=False, padding=0)
        self.box1.pack_start(self.button_convert, expand=True, fill=False, padding=0)
        self.box2.pack_start(self.text_input, expand=True, fill=True, padding=0)
        self.box2.pack_start(self.button_input, expand=True, fill=False, padding=0)
        self.box2.pack_start(self.check_input_open, expand=True, fill=False, padding=0)
        self.box3.pack_start(self.text_output, expand=True, fill=True, padding=0)
        self.box3.pack_start(self.button_output, expand=True, fill=False, padding=0)    
        self.box3.pack_start(self.check_output_open, expand=True, fill=False, padding=0)
        self.box4.pack_start(self.radio_dir, expand=True, fill=False, padding=0)    
        self.box4.pack_start(self.radio_file, expand=True, fill=False, padding=0)
        
    
    def show_all(self):
        self.combo_input.show()
        self.combo_output.show()
        self.button_convert.show()
        self.text_output.show()
        self.button_input.show()
        self.text_input.show()
        self.button_output.show()
        self.Vbox1.show()
        self.Hbox1.show()
        self.box1.show()
        self.separator.show()
        self.separator2.show()
        self.box2.show()
        self.box3.show()
        self.box4.show()
        self.check_input_open.show()
        self.check_output_open.show()
        self.radio_file.show()
        self.radio_dir.show()
        self.window.show()   
    ###################   Init Subfunctions  ####################################  


    def getActiveComboItem(self, combobox):
        model = combobox.get_model()
        active = combobox.get_active()
        if active < 0:
            return None
        return model[active][0]

    def convert(self, widget, data):
        inputFormat = self.getActiveComboItem(self.combo_input)
        outputFormat = self.getActiveComboItem(self.combo_output)
        convertor = formatConvertor.FormatConvertor()
        convertor.convert(self.text_input.get_text(), self.text_output.get_text(), inputFormat, outputFormat, 
                          "bla", 0, self.check_input_open.get_active(), self.check_output_open.get_active(), 0, 400000)

    def openDialog(self, widget, entry):
        if self.radio_dir.get_active():
            action = gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER
        else:
            action=gtk.FILE_CHOOSER_ACTION_OPEN
            
        self.chooser = gtk.FileChooserDialog(title=None,
                                             action=action,
                                             buttons=(gtk.STOCK_CANCEL,
                                                      gtk.RESPONSE_CANCEL,
                                                      gtk.STOCK_OPEN,
                                                      gtk.RESPONSE_OK))
        self.chooser.show()
        response = self.chooser.run()
        if response == gtk.RESPONSE_OK:
            print self.chooser.get_filename(), 'selected'
            if (action == gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER):
                entry.set_text(self.chooser.get_filename()+"/")
            else:
                entry.set_text(self.chooser.get_filename())
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
    
        
        self.chooser.destroy()
    
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False


    def dir_selected(self, obj, dir):
        """ The user has clicked on a directory on the left pane, so we need to load
        the files inside that dir on the right pane. """

        self.active_dir = dir
        self.populate_stop()

        self.stop_button.show()
        self.clear_button.set_sensitive(False)
        self.rename_button.set_sensitive(False)
        self.file_selected_model.clear()

        while gtk.events_pending():
            gtk.main_iteration()

        self.file_selected_model.clear()

        self.populate_selected_files(dir)
        self.selected_files.columns_autosize()
        self.statusbar.push(self.statusbar_context, _("Reading contents of dir %s") % dir)




def main():
    gtk.main()

if __name__ == "__main__":
    hello = PicosGtk()
    main()
