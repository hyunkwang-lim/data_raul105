Chapter 1: GRASP Contribution workflow {#chap01}
======================================

In this chapter you will get knowledge about the tools used to manage efficiently
the GRASP inversion code and the way to contribute the code. Here we talk about
how to send your contributions. In next chapter we will talk about types of contributions
and technical details to perform them.

The code is managed by GIT. It is a version control system widely-used. It
is open source and very powerful. If you don't know GIT, it is a good moment to start, 
once you know it you will never stop to use it. 

Additionally, GIT is accompanied by GitLab, a web interface to help users to access
the code and organize the community around it. GitLab offer a issue tracking, code reviewer 
and code contribution system (via pull request).

#### Index of content chapter 1:

- 1. [GRASP Contribution workflow](@ref chap01)
   - 1.1. [Code and repository organization]  (@ref chap0101)
   - 1.2. [Open an issue] (@ref chap0102)
   - 1.3. [Forking, contributing and pulling request] (@ref chap0103)

-------

# 1.1. Code and repository organization {#chap0101}

The GRASP code can be extended easily by external code. This external code is integrated
inside de GRASP code folder and compiled together. Following image represents the code organization
of GRASP where folders and repositories are mixed.

![Partial scheme of code folders and repositories](scheme-of-folder-and-repos.png)

This is also explained in user documentation: [www.grasp-open.com/doc/ch04.php#multi-repository](http://www.grasp-open.com/doc/ch04.php#multi-repository)

To help to prepare the code there is a tool called grasp-manager. The use of this tool
is explained in user documentation: [www.grasp-open.com/doc/ch04.php#grasp-manager](http://www.grasp-open.com/doc/ch04.php#grasp-manager)

From the technical point of view is important to identify all extensions in order to know where (in which repository)
the code is. Once that the repository is identified the user can contribute opening an issue or developing new code
via forking feature. Both things are explained in following sections.

# 1.2. Open an issue {#chap0102}

A way to support GRASP open is to follow the issue tracker. A user can help answering issues or
opening new one when a bug is identify. Creating an issue is an easy task, the difficult part is to provide
all required information to help the team to resolve the problem. Please, if you are going to open an 
issue try to provide as as much information as possible, such as a setting file and sdata file which help
to reproduce the problem. Last but not least, the user has to identify the proper repository where the issue 
has to be published. As it is explained in previous section, the GRASP code is organized with many repositories,
one for each extension. Please, before opening an issue identify the right repository. Following image 
shows the steps to open an issue:

![How to open an issue](issue-screenshot.png)

# 1.3. Forking, contributing and pulling request  {#chap0103}

The forking workflow is the way to contribute to GRASP open editing the repository.
The user can fix a bug, improve documentation or add new features. This is possible 
creating your own fork, modifying it and sending a pull request. This process is explained
in this guide of GitLab: [http://doc.gitlab.com/ee/workflow/forking_workflow.html](http://doc.gitlab.com/ee/workflow/forking_workflow.html)

