import pygame, sys, math


class Game:
    def __init__(self):
        self.FPS = 30
        self.windowWidth = 800
        self.windowHeight = 600

        # number of "blinks" for proximity cue
        self.proximityBlinks = 3
        
        self.isProximityCue = False
    
        pygame.init()
        self.display = pygame.display.set_mode((self.windowWidth, self.windowHeight))
        pygame.display.set_caption("Treasure Hunting in the Quantum Regime")
        self.clock = pygame.time.Clock()

        self.character = {image: None, rect: None, position: (0, 0), type: "Classical"}
        
        self.loadResources()

    def loadResources(self):
        #self.cueFilter = pygame.image.load(r'./res/cuebg.png')
        #self.cueFilterRect = self.cueFilter.get_rect()
        #self.cueFilter = pygame.transform.scale(self.cueFilter, (self.windowWidth, self.windowHeight))

        self.cueFilter = pygame.Surface((self.windowWidth,self.windowHeight), pygame.SRCALPHA)
        
        self.groundBg = pygame.image.load(r'./res/bg.png')
        self.groundBgRect = self.groundBg.get_rect()
        self.groundBg = pygame.transform.scale(self.groundBg, (self.windowWidth, self.windowHeight))

        self.character.image = pygame.image.load(r'./res/char.png')
        
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
                self.triggerProximityCue(0.0)

    def renderGraphics(self):
        self.display.blit(self.groundBg, (0, 0))

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
            self.renderGraphics()
            
            pygame.display.update()

            self.display.fill((0, 0, 0))
            self.clock.tick(self.FPS)

gameApp = Game()

if __name__ == "__main__":
    gameApp.run()
