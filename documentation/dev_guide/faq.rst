FAQ
---

* Is it ever ok to commit directly to master?

  Yes, but only under very restricted circumstances!  If in doubt make a branch
  and let someone else do the merge.

  + I've got a quick bugfix which I've tested - can I commit it to master?

    Well done on the testing.  A bugfix should go in a bugfix/XXX branch.  It's
    a single command to create this.  Another few commands and you'll have an pull
    request email to the hande-dev list for review.

  + But it's a really quick fix!  Surely it won't hurt?

    If it will affect functionality (and potentially someone else's jobs) then
    it probably ought to be reviewed!  If it's a very minor corner case of which
    you're certain, then commit to a bugfix branch and then do the merge
    yourself.  Always do this via a branch - don't commit directly to master.
    It's sensible to ask the original author if you're fixing their code
    however.

  + But I need to use this fix to make my runs work.

    You can always run from a bugfix branch.  Because you've committed it to the
    central git repository, you'll have access to it everywhere.

  + What if I need this bugfix to develop a new feature?

    I don't know.  Ask James!  One option is to base your subsequent feature
    branch off the bug fix branch before it's merged into master (git handles
    merges very well!) or to cherry-pick the bug fix into your feature branch or
    make enough noise to get the bug fix merged quickly.

  + I've added some comments to clear up something.

    This might be ok to commit to master.  If you designed the
    feature/documentation then you're effectively reviewing yourself.  If it's
    somebody else's code it's polite to have consulted someone on this (either
    by email, or a review branch).

  + But I've modified a feature that only I'm using...

    It sounds like this should be in an enhancement branch he/XXX.  If only
    you're using it it's even more important than someone else review it.

  + I've accidentally committed some changes to my local master.  What do I do?

    Remember that you can always push to a different branch on the main server.

    .. code-block:: bash

       $ git push origin master:he/XXX

    would push your changes to the he/XXX branch.  It's probably better, however
    to checkout your changes locally to a branch, and then roll back your
    master, and then commit the branch:

    .. code-block:: bash

       $ git checkout -b he/XXX
       $ git push --set-upstream origin he/XXX
       $ git checkout master
       $ git reset --hard origin/master

    Note the last command resets your local master to the same state as that on
    origin.  You should adapt the reset command to set your master to point to
    the desired commit (ie the first commit shared with the new branch he/XXX).

  + Ok - I've gone through the review process and I'd like to try to merge to
    master myself.  Is it easy?

    Easy as pie.  There's a workflow in the section Merging to master

* I've got a local branch which I've been working on for some time, but I don't
  want the pain of a large merge at the end.

  This sounds like a workflow problem.  Some comments on this:

  + We need to lose the idea of personal branches (note the branch namespace is
    organised by topic rather than person), even though a branch might be
    written entirely/mostly by one person.  In that sense, long-running
    development work should be split into small, logical chunks, each of which
    is attached one-at-a-time in its own branch.  We have always regretted
    having (multiple) long-running branches.
  + When wrenched away from a WIP with only a distant prospect of future free
    time, a commit and push with light notes is a very worthwhile thing.  It's
    probably even worthwhile committing a plan before committing any actual
    code.  If these are fast and flexible enough they will hopefully not
    discourage, but actually encourage organization.  It might also encourage
    (*gasp*) collaboration.  Perhaps you could create a directory in
    documentation as a place for such notes/roadmaps, somewhere between Python's
    PEP system and informal topic-based TODO lists?
  + We are pretty happy for development branches to be regularly rebased against
    master (*note*: not merged in either direction), to lessen the pain of one
    final merge between two very disparate branches.

* This is all very well (and I enjoy the Socratic method), but I'm stuck with
  a huge branch I don't have time to merge.  What do I do?

  Commit it as a feature/XXX or he/XXX and ask for help from the hande-dev list.

* How do I review code?

  We're working on a workflow for this.  One method is to make a branch (if
  you're not already in one) and just add comments to the source.  It's helpful
  if the review is part of the git history (even if the comments never actually
  make it to the master).  We currently are using `watson-style
  <http://goosecode.com/watson/>`_ tags in comments for code review and
  discussion, for example:

  .. code-block:: fortran

     ! [review] - JSS: How about doing it this way?
     ! [reply] - AJWT: I thought about it but that causes problems due to X.

  where JSS and AJWT are the initials of the reviewer and code author
  respectively.

* Will *my* code actually get reviewed?

  We're all usually terribly busy and have very little time, but in a group
  effort a little from each person goes a long way.  If you review others' code
  then they're more likely to review yours.  Make it easy to review, by keeping
  it clean and the features short.  Remember, this kind of review is far more
  lightweight than peer review of publications, and should be able to slot into
  people's 'free' time.  (Each branch is far more lightweight than a paper.)
  A simple pull-request should be enough to get people to review.  This is
  rather intricately tied in with the idea of project management.
  Prodding/cajoling/bullying emails are all possible to aid the review

* What happens if no-one replies to the pull request?

  Here are some opinions:

  + I suggest that after an agreed upon time (X working days?) without even
    a "I'll review but am too busy until next week" reply, the author is free to
    merge it into master (but should be open to fixes/improvements to that work
    that others subsequently suggest).
  + Having been burdened with years-long old dirty branches from other projects,
    merging is certainly vital.  I don't think lack of review should stop
    merging, but it should prompt someone to ask why.
  + I would view it as a sign that the work is stable and relatively
    complete (for the time being) and is ready to be used by others/in
    production calculations.

* What about major (long-term) development work?  Perhaps anyone engaged in
  major projects should send out 'pull-requests' to request review of ongoing
  work periodically?

  Yes.

* Why are we bothering with review?  Surely it makes life more difficult?

  In an attempt to avoid heaps of

  #. completely redundant code
  #. untested code
  #. buggy code

  all ending up in master.  The main reason is to encourage something resembling
  a coherent design and prevent someone going off in a (technical) direction
  others don't agree with/can see major problems with.  A big plus is that it
  helps everyone become familiar with code that they didn't write (which is why
  doing code review is good for newcomers).

* PhD students are going to be working on this. How do you see the work they
  produce on a single project over the course of 3 years going? How often should
  their code be subject to review?

  PhD projects are never one single monolithicproject (or at least shouldn't
  be!).  The amount and frequency of review is probably a function of how
  experienced a developer is (in general and with HANDE).  Remember a pull
  request can simply be an indication that the developer would like to start
  a conversation rather than presenting the final result.  Developers should
  also be encouraged to consider how a development task can be broken down into
  smaller projects, which might well aid design and testing, as well as reducing
  horrible merge conflicts from attempting to merge long-standing branches.

* How do I signify a 'fine - no need to comment' commit?

  We suggest a pull request to the email list followed immediately by an email
  announcing that the requester had also merged into master (or perhaps just the
  latter email).

