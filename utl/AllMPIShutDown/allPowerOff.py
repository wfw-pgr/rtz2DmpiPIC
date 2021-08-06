import subprocess
import getpass

def GetPassWord():
    passwd = getpass.getpass( 'PassWord :: ' )
    return passwd

def ClientPowerOff( host, passwd ):
    # -- SSH into Client and ps aux -- #
    cmd    = 'echo {0} | sudo -S shutdown -h now'.format( passwd )
    cmd    = 'ssh {0} {1}'.format( host, cmd )
    stdout = subprocess.call( cmd.split() )
    
if __name__ == '__main__':
    hosts  = [ 'kolmo02', 'kolmo03', 'kolmo04', 'kolmo05', \
               'kolmo06', 'kolmo07', 'kolmo08', 'kolmo09', 'kolmo10'  ]
    # hosts  = [ 'kolmo01', 'kolmo02', 'kolmo03', 'kolmo04', 'kolmo05', \
    #           'kolmo06', 'kolmo07', 'kolmo08', 'kolmo09', 'kolmo10'  ]
    passwd = GetPassWord()
    for host in hosts:
        ClientPowerOff( host, passwd )
