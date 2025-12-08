Current difficulties with AR:

1. There are weird things going on with points whose positions cannot be uniquely determined by the data we store.
   Consider the following examples:

- $\angle BOR=\angle YOD$; $B,D,O$ collinear; $\angle OBR=\angle BRO$; $\angle OBY=\angle BYO$, $\angle ORY=\angle RYO$.
  Where is $R$? It can either be on the same or opposite side of $BD$ with $Y$. This will lead to $\angle BOR=\angle BOY$. (I tackled this problem by both checking AR and **visual effects**: if they look like equal angles)
- $\angle BOR=\angle YOD$; $B,D,O$ collinear; $\angle OBR=\angle BRO$; $\angle OBY=\angle ORD$; $\angle OBY=\angle BYO$. Where is $D$? It can either be the same point as $B$ or not. In this case, if we are to add another relation $BR\perp DR$, then it will lead to $90\degree=0$. (I do not have a good idea for this).  
  Possible solution: such things happen because of division operations happening in checking the span in AR. For instance, $$2\angle A=2\angle B \implies \angle A=\angle B \text{ or } \angle B+90\degree \text{ or } \angle B -90\degree.$$
  The exact value should be derived from the diagram.

2. The current triangle congruency and similarity checks need to go through all triangles in every iteration, which significantly slows down the program.  
   Possible solution: the system can preprocess the diagram to get all possible pairs of congruent/similar triangles.
