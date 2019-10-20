import pygame, sys, math

WIN_W = 800
WIN_H = 600

class Game:
    def __init__(self, cellsize = [WIN_W/7, WIN_H/7]):
        self.FPS = 30

        self.cellsize = cellsize
        
        # number of "blinks" for proximity cue
        self.proximityBlinks = 3
        
        self.isProximityCue = False
    
        pygame.init()
        self.display = pygame.display.set_mode((WIN_W, WIN_H))
        pygame.display.set_caption("Treasure Hunting in the Quantum Regime")
        self.clock = pygame.time.Clock()

        self.character = Character(r'./res/char.png', "Classical", self)
        
        self.loadResources()
    def loadResources(self):
        #self.cueFilter = pygame.image.load(r'./res/cuebg.png')
        #self.cueFilterRect = self.cueFilter.get_rect()
        #self.cueFilter = pygame.transform.scale(self.cueFilter, (WIN_W, WIN_H))

        self.cueFilter = pygame.Surface((WIN_W,WIN_H), pygame.SRCALPHA)
        
        self.groundBg = pygame.image.load(r'./res/bg.png')
        self.groundBgRect = self.groundBg.get_rect()
        self.groundBg = pygame.transform.scale(self.groundBg, (WIN_W, WIN_H))

        self.character.loadResources(WIN_W//10, WIN_W//10)
        
    # proximity: float between [0, 1]
    def triggerProximityCue(self, proximity):
        self.isProximityCue = True
        self.currentProximityPhase = 2 * math.pi * self.proximityBlinks
        self.proximityCueColour = (proximity*255, 0, (1.0-proximity)*255)

    def handleEvents(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

            if event.type == pygame.MOUSEBUTTONDOWN:
                self.triggerProximityCue(1.0)

            if event.type == pygame.KEYDOWN:
                if not self.character.isMoving:
                    if event.key == pygame.K_LEFT:
                        self.character.triggerMovementCells("x", -1, 5)
                    if event.key == pygame.K_RIGHT:
                        self.character.triggerMovementCells("x", 1, 5)
                    if event.key == pygame.K_UP:
                        self.character.triggerMovementCells("y", -1, 5)
                    if event.key == pygame.K_DOWN:
                        self.character.triggerMovementCells("y", 1, 5)
                        
    def renderGraphics(self):
        # draw background
        self.display.blit(self.groundBg, (0, 0))

        # draw character
        self.display.blit(self.character.getImage(), self.character.position)
        if self.character.position[0]+self.cellsize[0] > WIN_W:
            self.display.blit(self.character.getImage(), ((self.character.position[0]+self.cellsize[0]) % WIN_W - self.cellsize[0], self.character.position[1]))

        if (self.character.position[1]+self.cellsize[1]) > WIN_H:
            self.display.blit(self.character.getImage(), (self.character.position[0], (self.character.position[1]+self.cellsize[1]) % WIN_H - self.cellsize[1]))

        if self.character.position[0] < 0:
            self.display.blit(self.character.getImage(), (self.character.position[0] + WIN_W, self.character.position[1]))

        if self.character.position[1] < 0:
            self.display.blit(self.character.getImage(), (self.character.position[0], self.character.position[1] + WIN_H))

        
        if self.isProximityCue:
            self.cueFilter.fill((*self.proximityCueColour, 255*(0.5 + 0.5*math.sin(self.currentProximityPhase))))
            self.display.blit(self.cueFilter, (0, 0))

            self.cueFilter.fill((*self.proximityCueColour, 255*(0.5 + 0.5*math.sin(self.currentProximityPhase))))
            
            #pygame.draw.rect(self.display,
            #                 (*self.proximityCueColour, 0.5 + 0.5*math.sin(self.currentProximityPhase)),
            #                 self.display.get_rect())
            self.currentProximityPhase -= math.pi/7.
            if self.currentProximityPhase <= 0:
                self.isProximityCue = False
            
    def run(self):
        while True:            
            self.handleEvents()
            self.character.update()
            self.renderGraphics()
            
            pygame.display.update()

            self.display.fill((0, 0, 0))
            self.clock.tick(self.FPS)

class Character:
    def __init__(self, imagedir, movementType, game, position = [0, 0]):
        self.image = None
        self.rect = None
        self.orientation = "right"
        self.position = position
        self.imagedir = imagedir
        self.movementType = movementType

        self.game = game
        
        self.isMoving = False

    def getCellPosition(self):
        return [(self.position[0]//self.game.cellsize[0]),
                (self.position[1]//self.game.cellsize[1])]

    
    def loadResources(self, w, h):
        self.image = pygame.transform.scale(pygame.image.load(self.imagedir), (w, h))
        self.images = {"left": pygame.transform.rotate(self.image, 180),
                       "right": self.image,
                       "up": pygame.transform.rotate(self.image, 90),
                       "down": pygame.transform.rotate(self.image, 270)}
        
        self.rect = self.image.get_rect()

    def triggerMovementToCell(self, axis, position1d, velocity=5):
        if axis == 'x':
            self.triggerMovement('x', (position1d - self.position[0]) * self.game.cellsize[0], velocity)
        else:
            self.triggerMovement('y', (position1d - self.position[1]) * self.game.cellsize[1], velocity)
            
    def triggerMovementCells(self, axis, cell_distance, velocity=5):
        if axis == 'x':
            self.triggerMovement('x', cell_distance * self.game.cellsize[0], velocity)
        else:
            self.triggerMovement('y', cell_distance * self.game.cellsize[1], velocity) 
            
    # axis is "x" or "y"
    def triggerMovement(self, axis, distance, velocity=5):
        self.isMoving = True
        self.distanceToMove = distance
        self.movementVelocity = velocity

        if axis == "x":
            self.orientation = ["left", "right"][int(distance > 0)]
            
        if axis == "y":
            self.orientation = ["up", "down"][int(distance > 0)]

    def update(self):
        if self.isMoving:
            if self.distanceToMove != 0:
                if self.orientation in ["down", "right"]:
                    self.position[[0, 1][{"right": 0, "down": 1}[self.orientation]]] += min(self.movementVelocity, self.distanceToMove)
                    self.distanceToMove = max(0, self.distanceToMove - abs(self.movementVelocity))
                    
                else:
                    self.position[[0, 1][{"left": 0, "up": 1}[self.orientation]]] += max(-self.movementVelocity, self.distanceToMove)
                    self.distanceToMove = min(0, self.distanceToMove + abs(self.movementVelocity))

                if self.distanceToMove == 0:
                    self.position[0] = self.position[0]%WIN_W
                    self.position[1] = self.position[1]%WIN_H
                    self.isMoving = False
                    
                    
    def getImage(self):
        return self.images[self.orientation]
            
gameApp = Game()

if __name__ == "__main__":
    gameApp.run()
