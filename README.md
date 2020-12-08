# HGTparadox
A repository containing the underlying code used in my 2021 HGT paper "The Gossip Paradox: why do bacteria share genes?"
The Paper will soon be released on the arXiv. This readme will be updated once suitable arXiv link exists, and again once the article is published in a pair reviewed journal.

"ManyCallScript_GrandSet.m" has the main guts of the agent based simulation. Agent based simulations use the magic of Binary tree information heaps in order to select locations and update probabilities quickly with each tick of step of our Markov chain- this does have the unfortunate side effect of making some of the code a little harder to read, but gives a **significant** speed up to simulations.

If anything is found to be missing, or if you have specific questions, please to not hesitate to contact me.
