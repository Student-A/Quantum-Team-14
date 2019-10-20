import pygame, sys, math

from buttons import *

WIN_W = 800
WIN_H = 600
GAME_W = 600

NUMBER_OF_CELLS = 9

class Game:
    def __init__(self, cellsize = [GAME_W/NUMBER_OF_CELLS, GAME_W/NUMBER_OF_CELLS]):
        self.gameState = {
            "p1_coins_left": 10, "p2_coins_left": 10,
            "treasure_cell_position": (2, 4)
        }
        
        self.FPS = 30

        self.cellsize = cellsize
        
        # number of "blinks" for proximity cue
        self.proximityBlinks = 3
        
        self.isProximityCue = False
        
        pygame.init()
        self.display = pygame.display.set_mode((WIN_W, WIN_H))
        self.gameDisplaySurf = pygame.Surface((GAME_W,GAME_W), pygame.SRCALPHA)

        # self.traceSurf = pygame.Surface((GAME_W,GAME_W), pygame.SRCALPHA)
        
        pygame.display.set_caption("Treasure Hunting in the Quantum Regime")
        self.clock = pygame.time.Clock()

        self.font = pygame.font.SysFont(None, 20, False, False)
            
        self.largeFont = pygame.font.SysFont(None, 24)

        self.hugeFont = pygame.font.SysFont(None, 40, False, False)

        pygame.mouse.set_cursor(*pygame.cursors.broken_x)

        self.coinFlipsLeftP1Text = Text("Coinflips left (CL): ", (610, 20), (50, 255, 50), self.largeFont)

        self.coinFlipsLeftP2Text = Text("Coinflips left (QM): ", (610, 60), (255, 50, 50), self.largeFont)

 

        
        self.coinFlipsP1CountText = Text(str(self.gameState["p1_coins_left"]), (770, 20), (255, 255, 255), self.largeFont)
        self.coinFlipsP2CountText = Text(str(self.gameState["p2_coins_left"]), (770, 60), (255, 255, 255), self.largeFont)
        
        self.character = Character(r'./res/char.png', "Classical", self, [self.cellsize[i] * (NUMBER_OF_CELLS-1)/2 + 1 for i in range(2)])

        self.updownText = Text("^ v", (620, 290), (255, 255, 255), self.hugeFont)
        self.leftrightText = Text("< >", (710, 290), (255, 255, 255), self.hugeFont)

        
        self.movementModeButton = RadioButton(10, (0, 255, 0), 2, 90, 0, 680, 300)

        self.plusCoinButton = AddButton((640, 380), (20, 20), (255, 255, 0))

        self.minusCoinButton = MinusButton((740, 380), (20, 20), (255, 255, 0))

        self.coinsToFlipText = Text("Coins to flip next move:", (610, 350), (255, 255, 255), self.largeFont)
        self.coinsWagerText = Text("0", (690, 380), (255, 0, 0), self.hugeFont)
        
        self.movementModeText = Text("Movement mode:", (610, 260), (255, 255, 255), self.largeFont)

        self.playButton = TextButton((670, 500), (0, 0, 255), "Play!", self.hugeFont)
        self.showProbabilityDistributionButton = TextButton((630, 570), (255, 0, 255), " Show Prob. Dist. ", self.largeFont)

        
        self.loadResources()

    def getProximityFromTreasure(self):
        distx = abs(self.gameState["treasure_cell_position"][0]-self.character.getCellPosition()[0])
        disty = abs(self.gameState["treasure_cell_position"][1]-self.character.getCellPosition()[1])

        distx = min(distx, NUMBER_OF_CELLS - distx)
        disty = min(disty, NUMBER_OF_CELLS - disty)

        return math.sqrt(distx**2 + disty**2)/(math.sqrt(0.5*NUMBER_OF_CELLS*NUMBER_OF_CELLS))
        
    def loadResources(self):
        #self.cueFilter = pygame.image.load(r'./res/cuebg.png')
        #self.cueFilterRect = self.cueFilter.get_rect()
        #self.cueFilter = pygame.transform.scale(self.cueFilter, (GAME_W, GAME_W))

        self.cueFilter1 = pygame.transform.scale(pygame.image.load(r'./res/cuebg.png'), (GAME_W, GAME_W))
        self.cueFilter2 = pygame.Surface((GAME_W,GAME_W), pygame.SRCALPHA)
        
        self.groundBg = pygame.image.load(r'./res/bg.png')
        self.groundBgRect = self.groundBg.get_rect()
        self.groundBg = pygame.transform.scale(self.groundBg, (GAME_W, GAME_W))

        self.character.loadResources(self.cellsize[0], self.cellsize[1])
        
    # proximity: float between [0, 1]
    def triggerProximityCue(self, proximity):
        self.isProximityCue = True
        self.currentProximityPhase = 2 * math.pi * self.proximityBlinks
        self.proximityCueColour = [proximity*255, 0, (1.0-proximity)*255]

    def handleEvents(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

            if event.type == pygame.MOUSEBUTTONDOWN:
                mousepos = pygame.mouse.get_pos()
                radioClicked = self.movementModeButton.checkButtonClicked(mousepos)
                
                if radioClicked:
                    pass
                
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

    def renderUI(self):
        self.coinFlipsLeftP1Text.render(self.display)
        self.coinFlipsLeftP2Text.render(self.display)

        self.coinFlipsP1CountText.render(self.display)
        self.coinFlipsP2CountText.render(self.display)

        self.movementModeButton.renderButtons(self.display)

        self.updownText.render(self.display)
        self.leftrightText.render(self.display)
        self.movementModeText.render(self.display)
        self.playButton.render(self.display)
        self.showProbabilityDistributionButton.render(self.display)


        
        self.plusCoinButton.render(self.display)
        self.minusCoinButton.render(self.display)
        self.coinsWagerText.render(self.display)

        self.coinsToFlipText.render(self.display)
        
    def renderGraphics(self):
        # draw background
        self.gameDisplaySurf.blit(self.groundBg, (0, 0))

        # character trace
        if len(self.character.movementHistory) >= 1:
            pygame.draw.line(self.gameDisplaySurf, (0, 255, 0),
                             [self.character.position[0]+self.cellsize[0]/2, self.character.position[1]+self.cellsize[1]/2],
                             [self.character.movementHistory[0][0]+self.cellsize[0]/2,
                              self.character.movementHistory[0][1]+self.cellsize[1]/2]
            )
        if len(self.character.movementHistory) >= 2:
            for i in range(len(self.character.movementHistory)-1):
                pygame.draw.line(self.gameDisplaySurf,
                                 (0, 150*(len(self.character.movementHistory)-i)/len(self.character.movementHistory)+50, 0),
                                [self.character.movementHistory[i][0]+self.cellsize[0]/2,
                                  self.character.movementHistory[i][1]+self.cellsize[1]/2],
                                 [self.character.movementHistory[i+1][0]+self.cellsize[0]/2,
                                  self.character.movementHistory[i+1][1]+self.cellsize[1]/2]
                )
        
        
        # draw character
        self.gameDisplaySurf.blit(self.character.getImage(), self.character.position)
        # character "wrap around" when moving near screen edges
        if self.character.position[0]+self.cellsize[0] > GAME_W:
            self.gameDisplaySurf.blit(self.character.getImage(), ((self.character.position[0]+self.cellsize[0]) % GAME_W - self.cellsize[0], self.character.position[1]))

        if (self.character.position[1]+self.cellsize[1]) > GAME_W:
            self.gameDisplaySurf.blit(self.character.getImage(), (self.character.position[0], (self.character.position[1]+self.cellsize[1]) % GAME_W - self.cellsize[1]))

        if self.character.position[0] < 0:
            self.gameDisplaySurf.blit(self.character.getImage(), (self.character.position[0] + GAME_W, self.character.position[1]))

        if self.character.position[1] < 0:
            self.gameDisplaySurf.blit(self.character.getImage(), (self.character.position[0], self.character.position[1] + GAME_W))

        
        # display proximity cue    
        if self.isProximityCue:
            self.cueFilter2.fill((*self.proximityCueColour, 255*(0.5 + 0.5*math.sin(self.currentProximityPhase))))

#            self.cueFilter1.fill((*self.proximityCueColour, 255*(0.5 + 0.5*math.sin(self.currentProximityPhase))))    

            self.gameDisplaySurf.blit(self.cueFilter2, (0, 0))

            
            self.currentProximityPhase -= math.pi/7.
            if self.currentProximityPhase <= 0:
                self.isProximityCue = False

        self.gameDisplaySurf.blit(self.cueFilter1, (0, 0))
            
        self.renderUI()
                
    def run(self):
        while True:            
            self.handleEvents()
            self.character.update()
            self.renderGraphics()

            self.display.blit(self.gameDisplaySurf, (0, 0))
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

        self.movementHistory = []
        self.movementHistorySize = 5
        
    def getCellPosition(self):
        return [(self.position[0]//self.game.cellsize[0]),
                (self.position[1]//self.game.cellsize[1])]
    
    def loadResources(self, w, h):
        self.image = pygame.transform.scale(pygame.image.load(self.imagedir), (int(w), int(h)))
        self.images = {"left": pygame.transform.rotate(self.image, 180),
                       "right": self.image,
                       "up": pygame.transform.rotate(self.image, 90),
                       "down": pygame.transform.rotate(self.image, 270)}
        
        self.rect = self.image.get_rect()

    def triggerMovementToCell(self, axis, position1d, velocity=5):
        if axis == 'x':
            self.triggerMovement('x', (position1d - self.getCellPosition()[0]) * self.game.cellsize[0], velocity)
        else:
            self.triggerMovement('y', (position1d - self.getCellPosition()[1]) * self.game.cellsize[1], velocity)
            
    def triggerMovementCells(self, axis, cell_distance, velocity=5):
        if axis == 'x':
            self.triggerMovement('x', cell_distance * self.game.cellsize[0], velocity)
        else:
            self.triggerMovement('y', cell_distance * self.game.cellsize[1], velocity) 
            
    # axis is "x" or "y"
    def triggerMovement(self, axis, distance, velocity=5):
        if self.isMoving:
            return
        
        self.isMoving = True
        self.distanceToMove = distance
        self.movementVelocity = velocity

   
        if len(self.movementHistory) > self.movementHistorySize:
            self.movementHistory = [self.position[:]] + self.movementHistory[:-1]
        else:
            self.movementHistory = [self.position[:]] + self.movementHistory[:]
            
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
                    self.position[0] = self.position[0]%GAME_W
                    self.position[1] = self.position[1]%GAME_W
                    self.isMoving = False
                    self.game.triggerProximityCue(1.0-self.game.getProximityFromTreasure())
                    
                    
    def getImage(self):
        return self.images[self.orientation]
            
gameApp = Game()

if __name__ == "__main__":
    gameApp.run()
