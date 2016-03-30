# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:40:11 2016

@author: cornkle
"""
import random

# Build a class
class Song(object):

# Build some modules for the class
# First thing is an initialisation module. Class knows "lyrics". 
    def __init__(self, lyrics):
        self.lyrics = lyrics

# What to do with the lyrics module? Print all its list parts. 
    def sing_me_a_song(self):
        for line in self.lyrics:
            print line
            
    def dance_for_me(self, randomword):
            self.lyrics.append(randomword)
            print sorted(self.lyrics, key=lambda *args: random.random())       

# Call class and feed it with lyrics.
happy_bday = Song(["Happy birthday to you",
                   "I don't want to get sued",
                   "So I'll stop right there"])
# Call class and feed it with lyrics. Again!
bulls_on_parade = Song(["They rally around tha family",
                        "With pockets full of shells"])
                        
just_something = Song(["This is an array", "Although in python a list", "I don't care"])                 

# Let your class do seomething. Let it do what you taught it
#happy_bday.sing_me_a_song()

#bulls_on_parade.sing_me_a_song()

#just_something.sing_me_a_song()

just_something.dance_for_me("fuckoff")