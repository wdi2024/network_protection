/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>

#include "header.h"  /* type definitions */
#include "input_graph.c" /* parser */
#include "timer.c"        /* timing routine */

#define MAX_LONG LONG_MAX

#define min(s, t) (s) < (t) ? (s) : (t)
#define max(s, t) (s) > (t) ? (s) : (t)







///////////////////////////////////////////
///////////////////////////////////////////The definition of functions
///////////////////////////////////////////


//The function for allocation
void *wpalloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}



RandomData* randinit(RandomData *rd )
{
  // long seed = (long)(timer() * 100000000);
  // printf("seed is %ld\n",seed);
  // srand(seed);

  rd->randNumIdx = 0;

  for (int i = 0; i < rd->maxLen; i++)
  {
    rd->randNums[i] = (cType)rand();
  }
  
  return rd;

}

RandomData* randinit_v(cType len)
{
  RandomData *rd = wpalloc(1, sizeof(RandomData));
  rd->maxLen = len;
  rd->randNums = (cType *)wpalloc(rd->maxLen, sizeof(cType));
  
  randinit(rd);

  return rd;
}

//The function for randomization
cType multi_rand(RandomData* rd)
{
  if(rd->randNumIdx >= rd->maxLen){
		randinit(rd);
  }	  
  return rd->randNums[rd->randNumIdx++];
}

/////////////////////////The two function for heap sorting edges 
void heapChange(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    // assert(idx[0] != idx[3]);
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void heapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        heapChange(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        // printf("swap 0 with %d \n",i);
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        heapChange(idx,edges,0,i-1);  
    }  

}  

#define CUT_SIGN 9
long number_of_cut_edge = 0;
  
/////////////////////The function to sort edges using capacity
void desendLinkWithRandomOrder(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    //pedges[i].tmp = -1*mrand(pd->rd) % (pedges[i].cap+1); //+1防止0的情况
    pedges[i].tmp = rand()%1000;
    if(pedges[i].w != CUT_SIGN){
      pedges[i].tmp += 10000; //非割边比割边大很多，因此排在最后，更倾向成为树割
    } 
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    heapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}

/////////////////////The function to sort edges using capacity
void desendLinkWithRandomOrderCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    //pedges[i].tmp = -1*mrand(pd->rd) % (pedges[i].cap+1); //+1防止0的情况
    pedges[i].tmp = 1000-pedges[i].cap+1; //+1防止0的情况
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    heapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}

///////////////////The function to sort edges using the value of currently minimal cut one edge belongs to 
void aesendLinkWithCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = multi_rand(pd->rd) % acv; 
    }
  }

    heapSort(idxs,pedges,cnt);
}

///////////////////按度数升序访问 
void aesendLinkWithDegree(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;
  nodeP* nodes = pd->gd->nodes;
  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    cType zn = pp->endNode;
    pedges[i].tmp =  -1 * (nodes+zn)->nIdx;//mrand(pd->rd) % acv; 
  }

    heapSort(idxs,pedges,cnt);
}

// //根据度数，但是同时
// void aOrderEdgeByDegreeAver(nodeP *np,PreprocData *pd)
// {
//   if (np->nIdx == 0)
//   {
//     return;
//   }
//   int cnt = np->nIdx;
//   nodeP* nodes = pd->gd->nodes;
//   sType *idxs = np->orderedEdges;
//   edgeP *pedges = np->edges;

//   for (int i = 0; i < cnt; i++)
//   {
//     edgeP *pp = pedges +i;
//     cType zn = pp->endNode;
//     pedges[i].tmp =  (nodes+zn)->nIdx;//mrand(pd->rd) % acv; 
//   }

//     HeapSort(idxs,pedges,cnt);
// }



/*
  the function to add more data to a traversal tree to accelerate the searching in the tree
  The idea is to precalcuate minimal cv value of a span of nodes in a traversal tree, e.g., when SPAN_LEN = 100, and a node has depth of 200, then the algorithm will pre-calculate the minimal cv value of the nodes between the node (dep=200) and an acestor(dep=101)
  upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
  lastDepMCV is the depth of the node depth that has the minimal cv in the span
  lastJointNodeId is the last ancestor node id that has more than one child nodes
  lastJointMCV is the cv of lastJoineNodeId

  改版后预处理算法要有变化：
      对每个节点的值，都要改进下考虑当前出发节点z的考虑下游分支的更小的cv'=cv2+-cof

  (1)处理时：
      在非段头节点中还要考虑父亲cv'的值，如果更小，则更新指向父亲
      跨段的cv2+-cof直接放到下段mcv中，因为solve中确认跨段了才使用当前段的mcv
  (2)求解时：
      看solve的注释


*/

long gRoot = 0;
void accComp(PreprocData *pd, cType curN, cType upid, cType upid_LC, long mcv, long* applyAdj, cType* depMap)
{
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  edgeP *pedges = (nodes + curN)->edges;

  long curCV = pd->gpcv[curN];
  cType curDep = pd->gpdep[curN];

  assert(curDep == 0 || curCV >0);

  depMap[curDep] = curN;
  pd->gpaccup[curN] = upid; //the previous segment tail

  int SPAN_LEN = pd->SPAN_LEN;
  int LN = pd->LEVEL_NUM;
  /******************************************/
  //作为祖先节点更新的点
  pd->gpaccup[curN] = upid; //the previous segment tail
  pd->gpaccposmcv[curN] = upid_LC; //the LN Child of previous segment tail
  
  if(curDep % SPAN_LEN == 0){
    upid = curN;
    upid_LC = 0; //节点id最小是1，所以0可以表示无效
  }
  else{
    //不需要更新什么，或者在后面逻辑中包含了

  }


  /****************************************/
  cType ances; 
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }

    ances = depMap[curDep - i];
    assert(ances != 0);
    if (applyAdj[ances] == 0) //ances的applyAdj是0，说明未被祖先设置过
    {
      if (pd->gcutadjsign[curN][i] == 1 ) // 
      {
        // printf("update ApplyAdj level i %ld, curN %ld dep %ld, to set applyAdj[%ld] dep %ld to %ld\n",i,curN,pd->gpdep[curN], ances,pd->gpdep[ances],pd->gcutadj[curN][i]);
        applyAdj[ances] = pd->gcutadj[curN][i];
        pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      }
    }

  }

  //更新完，此时curN的LN祖先(最上面祖先)值，在这个支线上肯定不会再变了，可以尝试更新minCV了
  if(curDep >= LN ){
      if(curDep % SPAN_LEN == LN){
        upid_LC = curN; //更新作为祖先的upid_LC值
        if(LN > 0){
          //如果LN是0，这就是上个段尾，不适用此更新规则
          pd->gpaccposmcv[curN] = upid_LC; //把当前节点LC也设置为最新，后面好处理 
        }
      }
      else if(curDep % SPAN_LEN == (LN+1)){
        mcv = MAX_LONG;
      }

      assert(applyAdj[ances] <= 0);
      assert(pd->gcutadj[ances][0] <= 0);
      assert(pd->gpcv[ances]+pd->gcutadj[ances][0] > 0);
      //例行更新祖先值，这个和段无关,而且是循环
      ances = depMap[curDep-LN];
      // printf("applyAdj[ances] %ld\n",applyAdj[ances]);
      //curN的poh装的是LN祖先的在包含curN及其小祖先的基础上的树的最小割，因为applyAdj是负值或0
      //需要反向补充回去的目标值，应该不是祖先的pcv,而是祖先pcv还要小的，树被裁剪后的最小值，应此时cv+[ances][0]
      pd->gpoh[curN] = pd->gpcv[ances] + pd->gcutadj[ances][0] - applyAdj[ances]; 
      mcv = min(mcv, pd->gpoh[curN]);
      pd->gpaccmcv[curN] = mcv; //这个值记录的是LN祖先的该段的从段首到该祖先的所有割值的最小

  }
  else{
    //nothing to update

  }


  while (cnt > 0)
  {
    cType zn = pedges->endNode;
    if (pd->gpfa[zn] == curN)
    {
      accComp(pd,zn, upid, upid_LC, mcv, applyAdj, depMap);
    }
    pedges++;
    cnt--;
  }

  //恢复添加到祖先身上的值，恢复成0
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }
    cType ances = depMap[curDep - i];
    // printf("curN %ld, dep %ld, ances %ld dep %ld\n",curN, pd->gpdep[curN], ances, pd->gpdep[ances]);
    // printf("level i %ld sign %ld, adj %ld  --  apladj %ld\n",i,pd->gcutadjsign[curN][i],pd->gcutadj[curN][i],applyAdj[ances] );
    assert(ances != 0);
    //这里是还结果的地方，之前赋值给applyAdj[ances]，这里还回去
    if (pd->gcutadjsign[curN][i] == 1 && pd->gcutadj[curN][i] > 0)
    {
      pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      assert(applyAdj[ances] == pd->gcutadj[curN][i]);
      applyAdj[ances] = 0;
    }
  }

}

void allCapSum(PreprocData *pd)
{
  nodeP *nodes = pd->gd->nodes;

  for (cType curN = 1; curN <= pd->gd->N; curN++)
  {
    nodeP *np = nodes + curN;
    edgeP *pedges = np->edges;
    int cnt = np->nIdx;
    np->totalCap = 0;
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + ni;
      np->totalCap += eh->cap;
    }
  }
}

void markPreProc(PreprocData *pd)
{
  // printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;

  for(cType curN =1; curN<=pd->gd->N; curN++){

      nodeP *np = nodes + curN;
      edgeP *pedges = np->edges;
      int cnt = np->nIdx;
      if (cnt < pd->K)
      {
        //少于K(不用等于)，这些都是割边
        for (int ni = 0; ni < cnt; ni++)
        {
          // nodeP* znp = nodes+eh->endNode;

          edgeP *eh = pedges + ni;
          if (eh->w != CUT_SIGN)
          {
            eh->rev->w = eh->w = CUT_SIGN; //标记为割边
            number_of_cut_edge++;
          }
        }
      }
  }
 

}


void stage1CutEdge(cType curN, PreprocData *pd)
{
  // printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  // printf("******curN is %ld\n",curN);
  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)wpalloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }
  sType *idxs = np->orderedEdges;

  // long cap;

  desendLinkWithRandomOrder(np,pd);

  // if (pd->mode == 1)
  // {
    // deOrderEdgeByRandomCap(np,pd); //
  // }
  // else if (pd->mode == 2)
  // {
  //   aOrderEdgeByAvgCV(np,pd);
  // }
  // else if (pd->mode == 3){
  //   //度数最小最小先访问
  //   aOrderEdgeByDegree(np,pd);
  // }

  // long sum_f = 0;

  
  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    // printf("zn is %ld (curN %ld) \n",zn,curN);
    if (zs == 1) // zn is an ancestor node or father node
    {
      // printf("\t zs is 1, zn %ld , *cruCV %ld += %ld\n",zn,*curCV,cap);
      *curCV += 1;
      pd->gpcv[zn] -= 1;

    }
    else if (zs == 0) // zn is not accessed, i.e., a child node
    {
      //poh临时使用，后面会作为cv2使用：ph记录curN的当前儿子zn子树遍历中，连到curN的边的容量的和，就是直接和curN连接的割值。所以需要置0
      // pd->gpoh[curN] = 0;

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      // printf("\t zs is 0, zn %ld , *cruCV %ld before markcut\n",zn, *curCV);
      stage1CutEdge(zn,pd);

      //assert(pd->gpoh[curN] > 0);
      // printf("----marCut return \n");
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      // printf("\t zs is %ld , zn %ld , *cruCV %ld += %ld after markcut\n",pd->gps[zn],zn, *curCV,pd->gpcv[zn]);
      *curCV += pd->gpcv[zn];

//       //此时：gpoh[curN]是zn子树连到curN的边总容量，即zn树所有节点脱离curN树后，curN树多出来的割
//         //zn树不连其他curN的子的树，只连curN或上面祖先
//       // pd->gpcv[zn] - pd->gpoh[curN]是zn子树在curN之上的割值, 即zn子树去掉后，curN减少的割值
		 // pd->gpoh[curN] 就是zn子树去掉后，增加的割值
         //cof_zn 就是zn树去掉后，整体增加的割值，如果小于0，就可以去掉
      // long cof_zn = pd->gpoh[curN] - ( pd->gpcv[zn] - pd->gpoh[curN]);
      // if(cof_zn < 0){ //如果去掉能进一步优化割值(减少)，则记录
      //   sum_f += cof_zn;
      // }
      // assert(pd->gpcv[zn] >= pd->gpoh[curN]);
      // pd->gpcof[zn] = cof_zn;
      // pd->gpcof[zn] = pd->gpoh[curN];

    }
    else
    {
      //这种情况就是,zn是curN的子树的叶子，正好连到curN有条边
      // printf("zs is %ld\n",zs);
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  // if (pd->gpoh[curN] == 0)
  // {
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n", cnt, curN, gRoot, *curCV, pd->gpdep[curN]);
  //   assert(1 == 2);
  // }
  //自己有没算进去？遍历所有邻居就是把自己算进去了
  //------------多层遍历算法步骤：每次遍历后，更新祖先增加的值到cut_adjust中
  

  // if(pd->gpoh[curN] == 0 && curN != gRoot){ //当邻居全是祖先是有可能的，已经在上面加了判断
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
  //   assert(1==2);
  // }
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }
  
  // assert(*curCV == 0 || *curCV > 1);
 

  // if(*curCV <= pd->K){
  //   printf("c CTY node %ld has cut value %ld\n",curN,*curCV);
  // }

  // assert(pd->gpoh[curN] >0);
    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  // printf("set curN %ld to stat 2 \n",curN);
  *curS = 2;

}


//
void stage2CutEdge(cType curN, PreprocData *pd, cType recent_depth )
{

// printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  if(*curCV < pd->K){
    recent_depth = *curDep;
  }

  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  sType *idxs = np->orderedEdges;

  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    if(pd->gpdep[zn] < recent_depth){
      //add the edge to set
      if(eh->w != CUT_SIGN){
        eh->rev->w = eh->w = CUT_SIGN; //标记为割边
        number_of_cut_edge ++;
      }
    }
    else if(pd->gpdep[zn] == *curDep + 1){
      //child
      stage2CutEdge(zn,pd,recent_depth);
    }
    

  }

}


/*

算法整体步骤：相对之前treem算法的改进
(1)预处理
	//按之前计算，加上一个值放到n上，代表对于n的father,如果去掉n这一支，f的mc变化
	sum_f = 0
	对f的每个为访问的child n:
		在遍历n前，在f上设置oh_f并置0
		在遍历计算中，每次访问祖先，除了计算mc，还更新oh_f (加上就可以)
		n返回后
			此时知道mc_n, oh_f(oh_f就是n这一支连到f的边的值，要包含f_n父子的边)
			这时计算n这一支去掉对f的mc的影响cof_n = oh_f - (mc_n-oh_f) = 2* oh_f - mc_n
				//oh_f是增加的割，mc_n-oh_f是之间到非f的割
			如果cof_n是负值 //说明割值降低了，值得庆贺
				就加到sum_f上,即sum_f 加上 负值，sum_f指的是f的所有减少割值的子n去掉，总共减少的割值
				如果不减少，这个n就不去掉
	//所有child遍历完之后计算mc2_f，就是f的最优的偏特树，此时这个树包含某些子树，而且mc是最小的
	得到mc2_f = mc_f + sum_f 
				
(2)计算时: solve 和 build的时候
	//因为mc2_f是已经去掉cof_n为负值的n的，正值就不管，还在mc_f中
	每次向上回溯，n回溯到f时，如果cof_n是负值，说明如果要保留n这一支，目前f最优割就包含n，即此时f用于计算的mc应该取(mc2_f + -cof_n)，即加上n这一支减掉的值
  现在问题来了：
    f的mc2和访问哪一支有关系，buildAcc预处理咋做？
    这样就意味着，不同的底层上来，每个节点的mc还不一样，导致 节点段 中最小值还不一样
    buildAcc记录的是向上的，所以可以记录的呀

*/
// long level_cumsum[110]; //保存子的每个level的负值的和
//////////////////////////////////The function to traverse the graph for one pass, i.e. the checkNode function of Algorithm 2 in the paper
void travForCut(cType curN, PreprocData *pd)
{
  // printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  // printf("******curN is %ld\n",curN);
  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)wpalloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedEdges;

  if (pd->mode == 1)
  {
    desendLinkWithRandomOrderCap(np,pd); //
  }
  else if (pd->mode == 2)
  {
    aesendLinkWithCV(np,pd);
  }
  else if (pd->mode == 3){
    //度数最小最小先访问
    aesendLinkWithDegree(np,pd);
  }

  // long sum_f = 0;

  //多层遍历算法步骤：每次遍历前，记录上面几位祖先的值到cut_adjust中
  //问题：还没有遍历成功
  cType fa = curN;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    if(pd->gpdep[fa] == 0){
      break;
    }
    // printf("fa is %ld, fasfa is %ld,  dep is %ld\n",fa,pd->gpfa[fa],pd->gpdep[fa]);
    fa = pd->gpfa[fa];
    // printf("   -> fa is %ld, dep is %ld\n",fa,pd->gpdep[fa]);
    
    pd->gcutadj[curN][i] = pd->gpoh[fa]; 
  }


  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    // printf("zn is %ld (curN %ld) \n",zn,curN);
    if (zs == 1) // zn is an ancestor node or father node
    {
      // printf("\t zs is 1, zn %ld , *cruCV %ld += %ld\n",zn,*curCV,cap);
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
      pd->gpoh[zn] += cap;
    }
    else if (zs == 0) // zn is not accessed, i.e., a child node
    {
      //poh临时使用，后面会作为cv2使用：ph记录curN的当前儿子zn子树遍历中，连到curN的边的容量的和，就是直接和curN连接的割值。所以需要置0
      // pd->gpoh[curN] = 0;

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      // printf("\t zs is 0, zn %ld , *cruCV %ld before markcut\n",zn, *curCV);
      travForCut(zn,pd);

      //assert(pd->gpoh[curN] > 0);
      // printf("----marCut return \n");
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      // printf("\t zs is %ld , zn %ld , *cruCV %ld += %ld after markcut\n",pd->gps[zn],zn, *curCV,pd->gpcv[zn]);
      *curCV += pd->gpcv[zn];

//       //此时：gpoh[curN]是zn子树连到curN的边总容量，即zn树所有节点脱离curN树后，curN树多出来的割
//         //zn树不连其他curN的子的树，只连curN或上面祖先
//       // pd->gpcv[zn] - pd->gpoh[curN]是zn子树在curN之上的割值, 即zn子树去掉后，curN减少的割值
		 // pd->gpoh[curN] 就是zn子树去掉后，增加的割值
         //cof_zn 就是zn树去掉后，整体增加的割值，如果小于0，就可以去掉
      // long cof_zn = pd->gpoh[curN] - ( pd->gpcv[zn] - pd->gpoh[curN]);
      // if(cof_zn < 0){ //如果去掉能进一步优化割值(减少)，则记录
      //   sum_f += cof_zn;
      // }
      // assert(pd->gpcv[zn] >= pd->gpoh[curN]);
      // pd->gpcof[zn] = cof_zn;
      // pd->gpcof[zn] = pd->gpoh[curN];

    }
    else
    {
      //这种情况就是,zn是curN的子树的叶子，正好连到curN有条边
      // printf("zs is %ld\n",zs);
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  // if (pd->gpoh[curN] == 0)
  // {
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n", cnt, curN, gRoot, *curCV, pd->gpdep[curN]);
  //   assert(1 == 2);
  // }
  //自己有没算进去？遍历所有邻居就是把自己算进去了
  //------------多层遍历算法步骤：每次遍历后，更新祖先增加的值到cut_adjust中
  fa = curN;
  cType faBelow = 0; //连fa到curN之前的容量值
  int actual_ln = -1;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    // printf("fa is %ld, dep is %ld\n",fa,pd->gpdep[fa]);
  
    if(pd->gpdep[fa] == 0){ //到root了，如果i=1说明自己就是root
	
      actual_ln = i-1;
      break;
    }

    fa = pd->gpfa[fa];
    //此时gcutadj保留的是curN中到fa这个祖先节点的边容量总和, 如果紧挨着，则包含上连边；否则不包含上连边
    pd->gcutadj[curN][i] = pd->gpoh[fa] - pd->gcutadj[curN][i];
    //(1) 下面计算后的faBelow是curN子树所有cut到这层祖先fa内部(包含fa)到curN前的(不包含curN)的边容量总和, faBelow包含curN的上连边
    faBelow += pd->gcutadj[curN][i];
    //(2) 而curN到fa之外(不包含fa)的cut值为  *curV - faBelow;
    //则：(1)-(2)的值是curN的子树全部去除后，对相应祖先fa的cut的增加值，负值更好; 
    pd->gcutadj[curN][i] = faBelow - (*curCV - faBelow);    
  
  }

  //actual_ln说明curN最大可以作为哪一层？ 0就是没有(curN就是root)
  if(actual_ln < 0){
    actual_ln = pd->LEVEL_NUM;
  }

  //如果level_num=5但是actual_ln=2,说明当前节点最多只能作为第2层，而且就是底层了
  //作为底层

  //作为最底层
  // if(pd->gcutadj[curN][actual_ln] < 0){
  //   pd->gcutadjsign[curN][actual_ln] = 1;
  // }

  //----------多层遍历算法步骤：

  /*
	下面记录curN对不同层的祖先，是否移除以及移除对祖先cv的改变
	如果actual_ln=1，说明curN的父亲就是root
  */
  // assert(actual_ln < 20);//给个限制，最多20层
  // memset(level_cumsum,0,(pd->LEVEL_NUM+1)*sizeof(long));
  //现在比较curN的i层和子的i+1层之间的值，两层的这两个值对应的祖先是同一个
  long tempCumSum = 0;
  pd->gcutadj[curN][0] = 0; //0层记录的是自己对自己的cv值的改变，或者说整个树中通过裁剪能达到的最大的可减少负值(保留自己为根的前提下)，该负值针对curN为根的树

  //特殊处理：对于L层祖先，要根据负值设置sign
  if(actual_ln == pd->LEVEL_NUM){
	  if(pd->gcutadj[curN][actual_ln] < 0){
		  pd->gcutadjsign[curN][actual_ln] = 1; //curN不会删除子孙了，如果<0就是移除他自己
	  }
  }
	
  //for (int i = actual_ln - 1; i >= 0; i--) //只能从curN的level-1开始，因为子是要level开始；i可以为0，因为对应子时i+1层; 
  for (int i = (actual_ln < pd->LEVEL_NUM ? actual_ln : actual_ln - 1); i >= 0; i--) //只能从curN的level-1开始，因为子是要level开始；i可以为0，因为对应子时i+1层; 如果actual_ln< L, 子可以比actual_ln大，否则得从L-1开始 
  {
	//对curN的i层祖先，计算子对该祖先的原始cv的负值的总和，存在这里
	//i=0时，就是对自己的影响
	tempCumSum = 0;

    //计算整个curN的树，若移除对curN的 i层祖先的CV的影响(就是子对其i+1层祖先)，所有负值保存在tempCumSum中；如果i=0，就是对自己的影响
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      assert(zn != 0);
      assert(zn != curN);

      if (pd->gpfa[zn] != curN)
      {
        continue;
      }

      if (pd->gcutadj[zn][i + 1] < 0)
      {
		  if(curN == 820812){ //TODEL including inside
			printf("curN %ld child zn %ld  val %ld sign %d \n",curN,zn, pd->gcutadj[zn][i+1],pd->gcutadjsign[zn][i+1]);
		  }
        tempCumSum += pd->gcutadj[zn][i + 1];
		
      }
    }

	if(curN == 820812){ //TODEL 包括里面
		printf("i %d, level_cumsum[i + 1] %ld\n",i,tempCumSum);
	}
	
    if (pd->gcutadj[curN][i] < tempCumSum) //tempCumSum是全部子的负值的总和(对祖先可减少的总负值)
    {
      // if(curN == 9904 && i == 4){
      //   printf("curN 9904 i=4 set sign self 1");
      // }
      pd->gcutadjsign[curN][i] = 1; //说明curN子树中，全树全部移除反而减少的更多
      // pd->gcutadjsign[curN][i+1] = 0; //设置的不是curN，应该是zn
    }
    else
    {


      // if(curN == 9904 && i == 4){
      //   printf("curN 9904 i=4 set sign children 1");
      // }
      pd->gcutadj[curN][i] = tempCumSum;
      pd->gcutadjsign[curN][i] = 0; //虽然用了子孙中所有负值，但是设0表示没有用自己全树参与
      // pd->gcutadjsign[curN][i+1] = 1;
    }

  }


  // if(pd->gpoh[curN] == 0 && curN != gRoot){ //当邻居全是祖先是有可能的，已经在上面加了判断
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
  //   assert(1==2);
  // }
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }

  // assert(pd->gpoh[curN] >0);

  if(pd->mode == 1){
//update ver 2 w,according to
    for (int ni = 0; ni < cnt; ni++)
    {
      // nodeP* znp = nodes+eh->endNode;

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      // nodeP *znp = nodes+zn;

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      //progate weight to curN's edges
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);//((eh->avgCV) * weight + *curCV)/(weight+1);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);//((reh->avgCV) * weight + *curCV)/(weight+1);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  // printf("set curN %ld to stat 2 \n",curN);
  *curS = 2;

}

////////////////////////////////The function to obtain min-cut value of given node pair, i.e., Algorithm 3 in the paper
/*
 (2)求解时：
      深度大的出发节点直接cv2，
      另外一个不是这个的祖先时，另外一个也可以用cv2
      [这个先不优化有点复杂]如果t是s祖先，
        如果cv2正好割掉下面，直接可以用cv2
        如果不是，也可以计算割掉对应树后这个t的割

      PS: t是s祖先，可以用cv2,t不是s的祖先也可以cv2，两个cv2都可以用
      问题是，如果t是s祖先，t的cv2对应的割并不能切断s这条支线，这样就不对了
        如果t是s祖先，其实我们算的是去掉这个支后的t的最小割
      那就简化：
        只要t不是s祖先，就可以用t的cv2

      后面的逐个逼近，需要用cv'
      
*/
//对给定节点，更新祖先值，并记录到按深度差为key的数组中，数组也就LN长度
void upByLevelNum(NodePropArr* pnp, cType startDep, cType curN, cType curDep, long *adj, cType LN)
{

  for (int i = 1; i <= LN; i++)
  {
    if ( curDep < (unsigned int)i || (curDep+LN) <= (startDep+i) )
    {
      break; //祖先已经不存在了，或已经到startDep的LN祖先了，没必要计算了
    }

    if (pnp->pacc_cut_adjust_sign[curN][i] == 1)
    {
      adj[startDep+i-(curDep)] = pnp->pacc_cut_adjust[curN][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
    }
  }
}

//函数只要检查，LN_v的最优值，是否包含v; 如果不包含，区分是v_top(不含)以下点最优不包含，还是v_top(含)以上
//
//v到LN_v的情况，每个点
long checkOptExclue(NodePropArr* pnp, cType v, cType v_top, cType LN_v, cType *pDep, cType *pFa, cType LN){
  
  cType minDep = 0;
  // printf("isOptimalExclude: v %ld dep %ld, LN_v %ld dep %ld\n",v,pDep[v],LN_v,pDep[LN_v]);
  while(v != LN_v){
    //如果有节点对应LN_v标志是1，说明肯定排除在最优之外，v的祖先有没有，v都被排除了。
    //如果v没有但祖先有，也可以通过再Fa回溯检测到
    if(pDep[v] - pDep[LN_v] <= LN && pnp->pacc_cut_adjust_sign[v][pDep[v] - pDep[LN_v]]==1){
      //如果不包含，还要体现到哪不包含
      minDep = pDep[v];    
    }
    v = pFa[v];
  }

  // printf("isOptimalExclude DONE \n");

  //除了v到v_top(不含)被包含或不包含的情况，还有一种是根本不知道，因为不在这个深度裁剪，这个时候什么情况？
  //这种咋处理?
  //从共同LN_v开始每个祖先，有几种可能：
  //(1)LN_v的最优包含一方，不包含另外一方
  //具体情况: 
  //(1) LN_s(不含)之下的祖先，单独看就行，其实就只有LN_s，其他也没意义
  /*(2) LN_s(含)及之上的祖先，则需要看不含一方的最高点，如果不含的最高点在LN_s(不含)以下，可以认为不含，因为去掉t不影响s。
        如果在之上就麻烦了，去掉t也去掉s了，这个就不好办了

  */
  //默认LN_v最优不需要去掉v
  //如果包含，那就是都包含，这个不用继续看了
  if(minDep > 0){
    if(minDep <= pDep[v_top]){
      return 2; //在v_top及之上断开了
    }
    else{
      return 1;
    }
  }

  return 0;
}

//其实还可以考察共同点v上面的祖先点W，只要v到w都没被丢弃，则s和t单独被丢弃的情况，都可以作为最优解来处理
//因为w不确定，其实是找最低的w, w的最小值
long calcuMinCutV4(PreprocData *pd, cType root, long minCandi, NodePropArr* pnp, cType s, cType t, int SPAN_LEN, aType LN)
{
  cType *pDep = pnp->pdep;
  long *pCV = pnp->pcv;
  cType *pFa = pnp->pfa;
  // long *poh = pnp->poh;
  // cType *paccup = pnp->pacc_upid;
  // cType *paccup_LC = pnp->pacc_pos_upmincv;
  // long *paccmcv = pnp->pacc_upmincv;

  // cType os = s, ot = t;
  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  printf("\nstart: s %ld dep %ld, t %ld dep %ld\n",s,pDep[s],t,pDep[t]);

  assert(pDep[s] >= pDep[t]);

  // long *adj = (long *)walloc(pd->gd->N+1, sizeof(long));
  long *adjs = (long *)wpalloc(pd->gd->N+1, sizeof(long));
  long *adjt = (long *)wpalloc(pd->gd->N+1, sizeof(long));
  // memset(adj, 0, (pd->gd->N+1) *  sizeof(long));
  memset(adjs, 0, (pd->gd->N+1) *  sizeof(long));
  memset(adjt, 0, (pd->gd->N+1) *  sizeof(long));

  long mcv = MAX_LONG;

  long startDeps = pDep[s];
  long startDept = pDep[t];
	
	//long temp;
  while(pDep[s] > pDep[t]){

    //先更新对祖先的补回值
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }
    }

    //调用祖先优化

	//temp = mcv; //TODEL
	
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
	
/*	
	if(mcv < temp){//TODEL 包括内部	
		printf("update using node s %ld , before mcv %ld\n", s,temp); 	
		printf("after mcv %ld\n", mcv); 
	}
*/
    s = pFa[s];

  }

  assert(pDep[s] == pDep[t]);

  if(s == t){
    goto end;
  }

  //此时s和t不相等，一起上溯
  while(s != t){
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }

      if (pnp->pacc_cut_adjust_sign[t][i] == 1)
      {
        adjs[ startDept - (pDep[t]-i) ] = pnp->pacc_cut_adjust[t][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }


    }

    //调用祖先优化
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
    mcv = min(mcv, pCV[t] + pnp->pacc_cut_adjust[t][0] - adjs[startDept - pDep[t]]); 

    s = pFa[s];
    t = pFa[t];
  }

end:
  free(adjs);
  free(adjt);
  return mcv;


}

//////////////////////function to load graph data
void inputGraph(PreprocData *pd){
  pd->gd = wpalloc(1,sizeof(GraphData));  

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes));

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", pd->gd->N, pd->gd->M);
}

///////////////////function to initialize data structure
void initializeAllData(PreprocData *pd){
  printf("c initPreProcData\n");
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)wpalloc(pd->total + 2, sizeof(cType));
  pd->allResults = wpalloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;

  allCapSum(pd);

  int LN = pd->LEVEL_NUM+1;
  for (int i = 0; i < pd->total; i++)
  {
    if(i>0){
      allResults[i] = allResults[0];
      continue;
    }

    allResults[i].pfa = (cType *)wpalloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)wpalloc(len, sizeof(cType));
    allResults[i].pcv = (long *)wpalloc(len, sizeof(long));
    allResults[i].poh = (long *)wpalloc(len, sizeof(long));
    allResults[i].pcof = (long *)wpalloc(len, sizeof(long));
    allResults[i].ps = (short *)wpalloc(len, sizeof(short));
    allResults[i].pacc_upid = (cType *)wpalloc(len, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)wpalloc(len, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)wpalloc(len, sizeof(cType));

    long* ptr = (long *)wpalloc(len*LN, sizeof(long));
    memset(ptr, 0, len * LN * sizeof(long));
    allResults[i].pacc_cut_adjust = (long **)wpalloc(len, sizeof(long*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust[j] = ptr+j*LN;
    }


    aType* ptr2 = (aType *)wpalloc(len*LN, sizeof(aType));
    memset(ptr2, 0, len * LN * sizeof(aType));
    allResults[i].pacc_cut_adjust_sign = (aType **)wpalloc(len, sizeof(aType*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust_sign[j] = ptr2+j*LN;
    }    

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].poh, 0, len * sizeof(long));
    memset(allResults[i].pcof, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    memset(allResults[i].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, len * sizeof(long));
    // memset(allResults[i].pacc_cut_adjust, 0, len * sizeof(cType*)); //存的指针不能初始化，指针再指向区域已经初始化了
    // memset(allResults[i].pacc_cut_adjust_sign, 0, len * sizeof(aType*));

  }  
}


/////////////////////function to traverse the graph data for multiple times, i.e., Algorithm 2 in the paper
void preProcessAll(PreprocData *pd){
  printf("c preProc \n");
  double tm;
  double totalProcTime = 0;
  NodePropArr *allResults = pd->allResults;
  //calculate total cap of one node
  cType root;

  cType len = pd->gd->N + 2;
  long *apply_adj = (long *)wpalloc(len, sizeof(long));
  cType *depth_map = (cType *)wpalloc(len, sizeof(cType));

  markPreProc(pd);
  printf("c total_cut_edge %ld after premark \n",number_of_cut_edge);

  for (int ipass = 0; ipass < pd->total; ipass++)
  {

    if(pd->rd != NULL){
      free(pd->rd->randNums);
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = randinit_v(pd->gd->M*2);    
    // printf("the %d times\n",i);
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gpoh = allResults[ipass].poh;
    pd->gpcof = allResults[ipass].pcof;
    pd->gps = allResults[ipass].ps;
    pd->gpaccup = allResults[ipass].pacc_upid;
    pd->gpaccmcv = allResults[ipass].pacc_upmincv;
    pd->gpaccposmcv = allResults[ipass].pacc_pos_upmincv;
    pd->gcutadj = allResults[ipass].pacc_cut_adjust;
    pd->gcutadjsign = allResults[ipass].pacc_cut_adjust_sign;

    if (pd->P == 300)
    {
      pd->mode = 3;
    }
    else
    {
      pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    }
	

    //----------------准备深度遍历   
    // if(pd->mode == 3){
    //   //此时预期是sf图，而根据我们生成sf图的方式，1是度数最大的节点
    //   root = 90000;
    // }
    // else{
    root = (rand() % pd->gd->N)+1;
    // }
    pd->roots[ipass] = root;
    pd->gpdep[root] = 0;
    printf("pass %d before markCut: root is %ld\n",ipass,root);
    fflush(stdout);
    //printf("pass %d, randidx %ld, root is %ld\n",i, randNumIdx,root);
    // printf("root fa %ld\n",allResults[i].pfa[root]);
    tm = timer();
    gRoot = root;
    
    stage1CutEdge(root,pd);
    stage2CutEdge(root,pd,0);

    pd->gpcv[root] = MAX_LONG;
    pd->gpoh[root] = MAX_LONG;


    //----------------准备构建加速结构
    /*
    printf("c before buildAcc\n");
	  fflush(stdout);

    memset(apply_adj, 0, len *  sizeof(long));
    memset(depth_map, 0, len *  sizeof(cType));
    buildAcc(pd, root, root, 0,MAX_LONG, apply_adj,depth_map);
    printf("c after buildAcc\n");
    fflush(stdout);
    */
    // printf("c 4 cv is %ld\n",pd->gpcv[4]);

    totalProcTime += timer() - tm;
    printf("c proctime for onepass: %10.06f\n", timer() - tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }

    memset(allResults[0].pfa, 0, len * sizeof(cType));
    memset(allResults[0].pdep, 0, len * sizeof(cType));
    memset(allResults[0].pcv, 0, len * sizeof(long));
    memset(allResults[0].poh, 0, len * sizeof(long));
    memset(allResults[0].pcof, 0, len * sizeof(long));
    memset(allResults[0].ps, 0, len * sizeof(short));
    memset(allResults[0].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[0].pacc_upmincv, 0, len * sizeof(long));

  }

  free(apply_adj);
  free(depth_map);

  printf("c preprocess times %10.6f\n", totalProcTime);
  printf("total_cut_edge %ld \n", number_of_cut_edge);

}

/////////////////////////The two function for heap sorting edges 
void heapChange2(sType *idx, edgeP ** edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    // assert(idx[0] != idx[3]);
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]]->tmp > edges[idx[i]]->tmp )    
            i++;  

        if(edges[idx[i]]->tmp <= edges[tempIdx]->tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void heapSort2(sType *idx, edgeP ** edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        heapChange2(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        // printf("swap 0 with %d \n",i);
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        heapChange2(idx,edges,0,i-1);  
    }  

} 

/////////////////////////////////////用来求解最小生成树的函数
////////只找非edge cut的边访问，即w == 9的边不能前进
edgeP** estack = NULL;
int stackTop = 0;
void spanningTreeBuild0(cType curN, PreprocData *pd, cType treeId)
{

// printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  pd->gpaccmcv[curN] = treeId;

  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + ni;
    cType zn = eh->endNode;
    if(eh->w == CUT_SIGN){ //割边
      //不访问另外一端，且只有两边都标记tree_id且不等时，再加入estack
      if(pd->gpaccmcv[zn] > 0 && pd->gpaccmcv[zn] != treeId){
        // printf("stack %ld is point %x \n",stackTop, eh);
        estack[stackTop++] = eh;
        // eh->tmp = eh->cap;
        eh->tmp = 1+(zn+curN)%10;
        // eh->tmp = (rand()%10)+1;
      }
    }
    else{
      if(pd->gpaccmcv[zn] == 0){
        //未访问过，访问之
        spanningTreeBuild0(zn,pd,treeId);
      }
      else{
        //已访问过，就不访问了
      }
    }

 
    

  }

}

sType* idxs = NULL;

//对所有找到的割边，结合其他边，通过最小生成树的方法找需要保护的边
void spanningTreeBuildAll(PreprocData *pd){
    //pd->gpaccmcv 记录所属的树
    //边的栈直接放边，可以找到双向边，就能得到所有顶点信息
    long M = pd->gd->M;
    long N = pd->gd->N;
    long len = M + 2;
  

    estack = (edgeP**)wpalloc(len, sizeof(edgeP*)); 
    memset(estack, 0, len * sizeof(edgeP*));
    stackTop = 0;

    idxs = (sType*)wpalloc(len,sizeof(sType));
    memset(idxs,0,len*sizeof(sType));  
    
    //buildMCT_cost0Tree 负责实现
    //深度遍历，躲开w=MarkSet的边，直到遍历返回
    //找未染色的点，继续上述过程，直到没有未染色的点
        //不一样的边，放入堆栈
    // cType treeId = 1;
    for(int i=1; i<=N; i++){
      if(pd->gpaccmcv[i] > 0){
        //已访问了，略过    
      }
      else{
        spanningTreeBuild0(i,pd, i);
      }
    }
    
    printf("stackTop is %d \n",stackTop);
    // printf("ggggg stack %ld is point %x \n",stackTop, estack[stackTop-1]);
    for(int i=0; i<stackTop; i++){
      idxs[i] = i;
    }

    // printf("ggggg2 stack %ld is point %x \n",stackTop, estack[stackTop-1]);
    //边从小到大排序
    heapSort2(idxs, estack, stackTop);  
    //For each 边 in stack 处理
    /*
      先根据映射树，确认两边节点所属于的最终的ID
        此时同时更新为id直接映射过去
      如果属于同个ID,则continue下一个
        否则可以合并树，得到仍然是树

        //合并意味着两个树id属于新id(也可以认为新set的集合),这个可以用映射表示(逻辑上近似于一棵树)

    */    
   long usedEdgeCount = 0;
   long totalCost = 0;
    for(int k=0; k<stackTop; k++){
      int i = idxs[k]; // 第k个位置就是排第k的estack中元素下标
      cType n1 = estack[i]->endNode;
      cType n2 = estack[i]->rev->endNode;
      cType t1 = n1;
      while(t1 != pd->gpaccmcv[t1]){
        assert(t1 > 0);
        t1 = pd->gpaccmcv[t1];
      }

      cType t2 = n2;
      while(t2 != pd->gpaccmcv[t2]){
        assert(t2 > 0);
        t2 = pd->gpaccmcv[t2];
      }

      pd->gpaccmcv[n1] = t1;
      pd->gpaccmcv[n2] = t2;

      if(t1 != t2){
        pd->gpaccmcv[t2] = t1;    
        usedEdgeCount ++;
        // printf("c use cut edge %ld %ld with cost %ld\n",n1, n2, estack[i]->tmp);
        totalCost += estack[i] ->tmp;      
      }
      else{
        //do nothing;
      }

    }

    printf("c total, chosen, cost is %ld\n", totalCost);

}





