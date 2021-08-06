#!/usr/bin/python
# -*- coding: utf-8 -*-
# Import smtplib for the actual sending function
import smtplib
import sys
import subprocess

# Import the email modules we'll need
from email.mime.text import MIMEText

# メイン関数
if __name__ == '__main__':

    # 以下の内容を変更する
    # me : 自分のGmail アドレス, you : 送信先のアドレス, passwd : Gmailパスワード
    job       = sys.argv[1]
    me        = "utplasmalec2017@gmail.com"
    passwd    = "MasaakiYamada"
    you       = "wfwknight@gmail.com"
    titletext = " [PIC] Nortification ---  job : " + job + " has just completed. [PIC] "
    body      = "PICシミュレーションが終了しました.以下，情報...\n"
    cmd       = "tail -n 150 job/" + job + "/run.log"
    tail      = subprocess.check_output( cmd.strip().split(" ") )
    body      = body + tail.decode("utf-8")

    msg = MIMEText(body)
    msg['Subject'] = titletext
    msg['From'] = me
    msg['To'] = you

    # Send the message via our own SMTP server.
    s = smtplib.SMTP('smtp.gmail.com',587)
    s.ehlo()
    s.starttls()
    s.ehlo()
    s.login(me, passwd)
    s.send_message(msg)
    s.close()
