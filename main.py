import pygame, sys

class Game:
    def __init__(self):
        self.FPS = 30
        
        pygame.init()
        self.display = pygame.display.set_mode((800, 600))
        pygame.display.set_caption("Treasure Hunting in the Quantum Regime")
        self.clock = pygame.time.Clock()

    def handleEvents(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
        
    def run(self):
        while True:
            self.handleEvents()

            pygame.display.update()
            self.clock.tick(self.FPS)

gameApp = Game()

if __name__ == "__main__":
    gameApp.run()
