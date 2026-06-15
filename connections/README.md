### **Server Access – EMBO POPGEN 2026**
Each student can access their assigned VM using their own username, hostname and password shared by slack. 
1) **Server Access via SSH** 
To connect to the VM, please open a terminal and use the following command:
```bash
ssh -i </path/to/private/key> <username>@<hostname>
````
For example, in the case of username "organizingsecretariat" and assigned hostname  popgen446-1.cloud.ba.infn.it, the command will be:
```bash
ssh -i </path/to/private/key> organizingsecretariat@popgen446-1.cloud.ba.infn.it
````
Please make sure that you use the private key corresponding to the public key previously shared for the setup of the VM access.

2) **Web platform access**
In addition to the SSH connection, students will also need to access two web platforms that will be used during the course.

They are available at the following addresses:

* `https://<hostname>`
and

* `https://<hostname>/jhub`
For example, if your assigned hostname is `popgen446-1.cloud.ba.infn.it`, the two links will be:

* https://popgen446-1.cloud.ba.infn.it
* https://popgen446-1.cloud.ba.infn.it/jhub
On both platforms, please log in using your personal username and password, shared via slack.

** Please keep your password strictly personal and do not share it with anyone **
