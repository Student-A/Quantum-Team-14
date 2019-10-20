import pygame
from pygame.locals import *
import math


def collide_w_mouse(x, y, recto):
    rect =pygame.Rect(recto)
    if (x > rect.left) and (x < rect.right) and (y > rect.top) and (y < rect.bottom):
        return True
    else:
        return False
    
def draw_button(window, bdy, buttontype, color, checked = False):
    if buttontype == "+":
        pygame.draw.rect(window, color, bdy, 2)
       
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))), 1)
        pygame.draw.line(window, (255,255,255), ((bdy[0]+(bdy[2]/2)), bdy[1]+2), ((bdy[0]+(bdy[2]/2)),bdy[1]+bdy[3]-2), 1)
    if buttontype == 1:
        pygame.draw.rect(window, color, bdy, 2)
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))))
    if buttontype == "++":
        pygame.draw.rect(window, color, bdy, 2)
       
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+(bdy[2]/2)-2,(bdy[1]+(bdy[3]/2))), 1)
        pygame.draw.line(window, (255,255,255), ((bdy[0]+(bdy[2]/4)), bdy[1]+2), ((bdy[0]+(bdy[2]/4)),bdy[1]+bdy[3]-1), 1)
 
        pygame.draw.line(window, (255,255,255), (bdy[0]+(bdy[2]/2)+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))), 1)
        pygame.draw.line(window, (255,255,255), ((bdy[0]+((bdy[2]/2)+(bdy[2]/4))), bdy[1]+2), ((bdy[0]+((bdy[2]/2)+(bdy[2]/4))),bdy[1]+bdy[3]-1), 1)
 
    if buttontype == "--":
        pygame.draw.rect(window, color, bdy, 2)
       
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+(bdy[2]/2)-2,(bdy[1]+(bdy[3]/2))), 1)
        pygame.draw.line(window, (255,255,255), (bdy[0]+(bdy[2]/2)+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))), 1)

    if buttontype == "checkbox":
        pygame.draw.rect(window, color, bdy, 2)
        if checked:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)
    """if buttontype == 5:
        pygame.draw.rect(window, color, bdy, 2)
        if DEFLECT == True:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)
    if buttontype == 6:
        pygame.draw.rect(window, color, bdy, 2)
        if LONGITUDINAL == True:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)
    if buttontype == 7:
        pygame.draw.rect(window, color, bdy, 2)
        pygame.draw.line(window, (255,255,255), (bdy[0]+8, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-8,(bdy[1]+(bdy[3]/2))), 6)
    if buttontype == 8:
        pygame.draw.rect(window, color, bdy, 2)
        pygame.draw.line(window, (255,0,0), (bdy[0]+5, bdy[1]+2), (bdy[0]+bdy[2]-5,(bdy[1]+bdy[3]-2)), 2)
        pygame.draw.line(window, (255,0,0), ((bdy[0]+5), bdy[1]+bdy[3]-2), (bdy[0]+bdy[2]-5, bdy[1]+2), 2)
    if buttontype == 9:
        pygame.draw.rect(window, color, bdy, 2)
        if SHOW_B == True:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)
    if buttontype == 10:
        pygame.draw.rect(window, color, bdy, 2)
        if PAUSE == False:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)
    """

class RadioButton:
    def __init__(self, size, color, choices, x_distance, y_distance, x, y):
        self._size=size
        self._choices=choices
        self._x_distance=x_distance
        self._y_distance=y_distance
        self._color=color
        self._current_choice=0
        self._position_x=x
        self._position_y=y
        

    def renderButtons( self, surface):
        x = self._position_x
        y = self._position_y
        
        for button_count in range(self._choices):
            pygame.draw.circle(surface, self._color, (x+(self._x_distance*button_count), y+(self._y_distance*button_count)), self._size, 2)
            if (button_count==self._current_choice):
                if (self._size-3 > 0):
                    highlight_size=self._size-3
                else:
                    highlight_size=1
                pygame.draw.circle(surface, (255,0,0), (x+(self._x_distance*button_count), y+(self._y_distance*button_count)), highlight_size,4)

    def checkButtonClicked( self, mouse_position ):
        for button_count in range(self._choices):
            if (math.sqrt((mouse_position[0]-(self._position_x+self._x_distance*button_count))**2 + (mouse_position[1]-(self._position_y+self._y_distance*button_count))**2) <= self._size):
                self._current_choice=button_count
                return str(button_count) 
        return False
        

class Button:
    
    def __init__(self, position, size, colour):
        self.rect = pygame.Rect(position[0], position[1], size[0], size[1])
        self.position = position
        self.size = size
        self.colour = colour
        
    def isPointInside(self, x, y):
        return (x > self.rect.left) and (x < self.rect.right) and (y > self.rect.top) and (y < self.rect.bottom)

    def render(self):
        pass
    
class AddButton(Button):
    def __init__(self, position, size, colour):
        Button.__init__(self, position, size, colour)
        
    def render(self, window):
        pygame.draw.rect(window, self.colour, self.rect, 2)

        bdy = list(self.rect)
        
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))), 1)
        pygame.draw.line(window, (255,255,255), ((bdy[0]+(bdy[2]/2)), bdy[1]+2), ((bdy[0]+(bdy[2]/2)),bdy[1]+bdy[3]-2), 1)
        

class MinusButton(AddButton):        
    def render(self, window):
        bdy = list(self.rect)
        
        pygame.draw.rect(window, self.colour, bdy, 2)
        
        pygame.draw.line(window, (255,255,255), (bdy[0]+2, (bdy[1]+(bdy[3]/2))), (bdy[0]+bdy[2]-2,(bdy[1]+(bdy[3]/2))))    
    
class CheckboxButton(Button):
    def __init__(self, position, size, colour, checked = False):
        Button.__init__(self, position, size, colour)
        self.isChecked = checked
        
    def render(self, window):
        bdy = list(self.rect)
        pygame.draw.rect(window, self.colour, bdy, 2)
        if self.isChecked:
            pygame.draw.lines(window, (255,255,255), False, [(bdy[0]+2, bdy[1]+2), (bdy[0]+6, bdy[1]+bdy[3])  ,(bdy[0]+bdy[2], bdy[1])], 2)    

    def toggle(self):
        isChecked = not isChecked

class TextButton(Button):
    def __init__(self, position, colour, text, font):
        
        self.text = text
        self.textObj = Text(text, position, colour, font)

        self.rect = self.textObj.rect
        self.rect.x = position[0]
        self.rect.y = position[1]

        Button.__init__(self, position, (self.rect.w, self.rect.h), colour)
        
    def render(self, window):
        pygame.draw.rect(window, self.colour, self.rect, 2)
        self.textObj.render(window)
        
        
class Text:
    def __init__(self, text, position, colour, font):
        self.text = text
        self.position = position
        self.colour = colour
        
        self.surf = font.render(text, True, colour)
        self.rect = self.surf.get_rect()
        self.rect.x = position[0]
        self.rect.y = position[1]

    def render(self, display):
        display.blit(self.surf, self.rect)
        
    def changeText(self, newText, newColour):
        self.text = newText
        self.colour = newColour
        self.redraw()

    def redraw(self):        
        self.surf = font.render(self.text, True, self.colour)
        self.rect = self.surf.get_rect()
        self.rect.x = self.position[0]
        self.rect.y = self.position[1]

