Current difficulties with AR:

1. There are weird things going on with points whose positions cannot be uniquely determined by the data we store. 
Consider the following examples:

- $\angle BOR=\angle YOD$; $B,D,O$ collinear; $\angle OBR=\angle BRO$; $\angle OBY=\angle BYO$, $\angle ORY=\angle RYO$. 
Where is $R$? It can either be on the same or opposite side of $BD$ with $Y$. This will lead to $\angle BOR=\angle BOY$. (I tackled this problem by both checking AR and __visual effects__: if they look like equal angles)

- $\angle BOR=\angle YOD$; $B,D,O$ collinear; $\angle OBR=\angle BRO$; $\angle OBY=\angle ORD$; $\angle OBY=\angle BYO$. Where is $D$? It can either be the same point as $B$ or not. In this case, if we are to add another relation $BR\perp DR$, then it will lead to $90\degree=0$. (I do not have a good idea for this)

2. The traceback function of AR needs further investigation. Specifically, the order of adding parents to the relations list, checking duplicates, and adding the child relation matters. 