#include "bolt/Passes/PostLinkOutliningPass.h"
#include "bolt/Core/BinaryContext.h"
#include "bolt/Passes/BinaryPasses.h"
#include "bolt/Core/BinaryFunction.h"
#include "bolt/Passes/DataflowInfoManager.h"
#include "bolt/Core/MCPlus.h"
#include "llvm/Support/CommandLine.h"
#include "bolt/Core/HashUtilities.h"
#include "MCTargetDesc/AArch64MCTargetDesc.h"

namespace llvm {
namespace bolt {

bool PostLinkOutlining::isCalleeAccessingCallerStack(const MCInst &Inst,BinaryFunction *Callee, BinaryContext &BC) {
  /*if (!Callee || Callee->empty()) return false;
    const auto &Desc = BC.MII->get(Inst.getOpcode());
    bool IsLoad = Desc.mayLoad();
  // 扫描 Callee 的入口基本块
    if (Callee->getLayout().block_empty()) return false;
    BinaryBasicBlock &EntryBB = *Callee->getLayout().block_front();
  for (auto &Inst : EntryBB) {
    // 如果发现了修改 SP 的指令（通常是 sub sp, sp, #imm），说明进入了它自己的栈帧管理
    //if (BC.MIB->isStackPointerModification(Inst)) break;

    // 如果在修改 SP 之前，发现了读取 SP 正偏移的操作
    // 例如：ldr x0, [sp, #0] 或 ldr x1, [sp, #8]
    bool IsStore, IsPrologue, IsEpilogue;
uint16_t Reg;
int64_t Offset;
uint8_t Size;
bool IsStp, IsLdp;
MCPhysReg BaseReg;

// 调用 11 个参数的版本
if (IsLoad&&BC.MIB->isStackAccess(Inst, IsLoad, IsStore, IsPrologue, IsEpilogue, 
                          BaseReg, Reg, Offset, Size, IsStp, IsLdp)) {
      int64_t Offset;
      // 这里的逻辑需要检查操作数，看偏移量是否为正
      // 只要在它自己开辟空间前读取了 sp，就说明它在向调用者要参数
      return true; 
    }
  }
  return false;*/
}
bool PostLinkOutlining::isInShrinkWrappedRegion(const BinaryBasicBlock &BB,const MCInst &Inst, const BinaryFunction &BF) {
  /*
  // 1. 获取该指令所属的基本块 (BinaryBasicBlock)
  const BinaryBasicBlock *BB = BF.getBasicBlockForInst(&Inst);
  
  // 如果指令不属于任何基本块（通常不应该发生），返回 false
  if (!BB) {
    return false;
  }

  // 2. 检查该基本块是否在收缩包装的范围内
  // 在 BOLT 中，基本块通常包含一个状态位或属性来标识它是否被“包装”
  // 即：该块执行时，新的 Prologue 已经执行，而 Epilogue 尚未执行
  return BB->isShrinkWrapped();*/
}
bool PostLinkOutlining::areSequencesEquivalent(const std::vector<MCInst> &Seq1, 
                                               const std::vector<MCInst> &Seq2,
                                               const BinaryContext &BC) {
  if (Seq1.size() != Seq2.size()) return false;

  // 定义符号对比逻辑：在 Outlining 中，通常要求符号引用必须完全相同
  auto SymbolComp = [](const MCSymbol *S1, const MCSymbol *S2) {
    return S1 == S2;
  };

  for (size_t i = 0; i < Seq1.size(); ++i) {
    // 直接使用 BC.MIB->equals
    if (!BC.MIB->equals(Seq1[i], Seq2[i], SymbolComp)) {
      return false;
    }
  }
  return true;
}
bool PostLinkOutlining::hasMoreThanEightArgs(BinaryBasicBlock &BB, 
                                            BinaryBasicBlock::iterator CallIt, 
                                            BinaryContext &BC) {
  MCInst &CallInst = *CallIt;
  BinaryFunction *Callee = BC.getFunctionForSymbol(BC.MIB->getTargetSymbol(CallInst));

  // 维度一：向内看（分析目标函数）
  if (isCalleeAccessingCallerStack(CallInst,Callee, BC)) return true;

  // 维度二：向回看（分析调用前的指令）
  auto It = CallIt;
  size_t Window = 10; // 往前看10条指令
  while (It != BB.begin() && Window--) {
    --It;
    // 如果发现 str ..., [sp, #...] 且偏移 >= 0
   /* if (BC.MIB->isStore(*It) && BC.MIB->isStackAccess(*It)) {
      // 识别到栈传参准备动作
      return true;
    }*/
    // 如果碰到了另一个 call，说明进入了另一个调用的范畴，停止回溯
    if (BC.MIB->isCall(*It)) break;
  }

  return false;
}
std::vector<PostLinkOutlining::SequenceInstance> PostLinkOutlining::getAllseqs(BinaryContext &BC, uint64_t len) {
    std::vector<PostLinkOutlining::SequenceInstance> AllSeqs;

    for (auto &BFI : BC.getBinaryFunctions()) {
        BinaryFunction &BF = BFI.second;
        //if (!shouldOptimize(BF)) continue;
        bool isLeaf=true;
        for (const BinaryBasicBlock &BB : BF) {
          for (const MCInst &Inst : BB) {
          // 检查是否为 Call 指令
            if (BC.MIB->isCall(Inst)) {
              isLeaf=false; 
            }
          }
        }
        if (!BF.hasConstantIsland() && isLeaf) {
          continue;
        }
        for (BinaryBasicBlock &BB : BF) {
            if (BB.size() < len) continue;

            // 滑动窗口遍历基本块内的所有指令
            size_t CurrentIdx = 0;
            for (auto I = BB.begin(), E = BB.end(); std::distance(I, E) >= (int64_t)len; ++I, ++CurrentIdx) {
                
                std::vector<MCInst> CurrentInsts;
                auto WindowIt = I;
                for (uint64_t k = 0; k < len; ++k, ++WindowIt) {
                    CurrentInsts.push_back(*WindowIt);
                }

                // 检查序列是否可以被安全外联（无分支指令等）
                if (isOutlinable(BF,BB,CurrentInsts, BC)) {
                    SequenceInstance Instance;
                    Instance.Insts = std::move(CurrentInsts);
                    //Instance.Hash = computeHash(Instance.Insts, BC);
                    Instance.BF = &BF;
                    Instance.BB = &BB;
                    Instance.It = I;
                    Instance.StartIndex = CurrentIdx;

                    AllSeqs.push_back(Instance);
                    /*std::advance(I, len-1); 
                    CurrentIdx += len-1;*/
                }
            }
        }
    }
    return AllSeqs;
}
void PostLinkOutlining::setLabel(SequenceInstance &Seq) {
    // 在内存记录中将其标记为已使用
    Seq.IsLabeled = true;
    
}
BinaryFunction* PostLinkOutlining::createOutlinedFunction(
    BinaryContext &BC, 
    const std::vector<MCInst> &Seq, 
    unsigned Index){
  std::string FuncName = "outline_" + std::to_string(Index);
  BinaryFunction *OutlinedFunc = BC.createInjectedBinaryFunction(FuncName, true);
  std::vector<MCInst> Instrs;
  const auto *MIB = BC.MIB.get();
  const MCPhysReg FP = MIB->getFramePointer();
  const MCPhysReg LR = BC.MRI->getRARegister();
  MCInst StpInst;
  StpInst.setOpcode(AArch64::STPXpre);
  StpInst.addOperand(MCOperand::createReg(AArch64::SP));      // 预索引结果写回 SP
  StpInst.addOperand(MCOperand::createReg(FP));      // X29
  StpInst.addOperand(MCOperand::createReg(LR));      // X30
  StpInst.addOperand(MCOperand::createReg(AArch64::SP));      // 基址
  StpInst.addOperand(MCOperand::createImm(-2));      // 偏移: -2 * 8 = -16 bytes
  Instrs.push_back(StpInst);
  for (const auto &Inst : Seq) {
    Instrs.push_back(Inst);
  }
  MCInst LdpInst;
  LdpInst.setOpcode(AArch64::LDPXpost);
  LdpInst.addOperand(MCOperand::createReg(AArch64::SP));      // 后索引结果写回 SP
  LdpInst.addOperand(MCOperand::createReg(FP));      // X29
  LdpInst.addOperand(MCOperand::createReg(LR));      // X30
  LdpInst.addOperand(MCOperand::createReg(AArch64::SP));      // 基址
  LdpInst.addOperand(MCOperand::createImm(2));       // 偏移: 2 * 8 = 16 bytes
  Instrs.push_back(LdpInst);
// 通用返回指令
  MCInst RetInst;
  MIB->createReturn(RetInst);
  Instrs.push_back(RetInst);

  std::vector<std::unique_ptr<BinaryBasicBlock>> BBs;

  BBs.emplace_back(OutlinedFunc->createBasicBlock());
  BBs.back()->addInstructions(Instrs.begin(), Instrs.end());
  BBs.back()->setCFIState(0);
  OutlinedFunc->insertBasicBlocks(nullptr, std::move(BBs), true, false);
  OutlinedFunc->setAlignment(BC.isAArch64() ? 4 : 16);
  OutlinedFunc->updateState(BinaryFunction::State::CFG_Finalized);

  return OutlinedFunc;  
}
/*bool PostLinkOutlining::isOutlinable(BinaryFunction &BF,BinaryBasicBlock &BB,const std::vector<MCInst> &Seq, BinaryContext &BC) {
for (const MCInst &Inst : Seq) {
        const MCInstrDesc &Desc = BC.MII->get(Inst.getOpcode());

        // 1. 终止符与返回指令检查
        // 序列内部绝不能包含跳转、分支或返回，必须保持控制流线性
    if ((Desc.isTerminator() || Desc.isReturn() || BC.MIB->isBranch(Inst)) &&
        !BC.MIB->isCall(Inst)) {
        return false;
    }
    // 2. 检查 LR (X30) 的使用和定义 (手动遍历操作数以确保兼容性)
    if (!BC.MIB->isCall(Inst)) {
        // 检查所有寄存器操作数，涵盖 Def 和 Use
        for (unsigned i = 0; i < Inst.getNumOperands(); ++i) {
            const MCOperand &Op = Inst.getOperand(i);
            if (Op.isReg() && Op.getReg() == llvm::AArch64::X30) {
                return false; // 只要碰了 X30 且不是 Call，就不合格
            }
        }   
    
        // 额外检查隐式定义的寄存器 (针对某些特殊的系统指令)
        const MCInstrDesc &Desc = BC.MII->get(Inst.getOpcode());
        if (Desc.hasImplicitDefOfPhysReg(llvm::AArch64::X30) ||
            Desc.hasImplicitUseOfPhysReg(llvm::AArch64::X30)) {
            return false;
        }
    }
        // 3. 论文核心逻辑：封装函数调用 (Wrapped Function Calls)
        // 论文提到：如果包含 bl 指令，必须检查其参数传递方式
        if (BC.MIB->isCall(Inst)) {
            // AArch64 规定前 8 个参数通过 x0-x7 传递
            // 如果调用涉及栈上传参（超过 8 个参数），会修改当前函数的栈帧
            // 这会导致外联后的栈偏移计算极其复杂，因此视为不合格
            if (hasMoreThanEightArgs(BB,&Inst, BC)) {
                return false;
            }
            
            // 另外，在 Outline 序列中间有 Call 会破坏 LR (x30)
            // 我们的处理逻辑需确保在新函数内正确保存/恢复 x30
        }



        // 4. PC 相关指令检查 (AArch64 ADR/ADRP)
        // 因为外联后物理地址改变，所有的 PC-relative 计算都会失效
        if (BC.MIB->isPCRelative(Inst)) {
            return false;
        }

        // 5. 论文避让逻辑：收缩包装 (Shrink-wrapping)
        // 如果该指令位于已被编译器进行收缩包装优化的区域，外联会破坏其优化效果
        // 在 BOLT 中可通过检查基本块是否在 Prologue/Epilogue 移动路径上来识别
        if (isInShrinkWrappedRegion(BB,Inst, BF)) {
            return false;
        }

        // 6. 伪指令与系统指令排除
        if (BC.MIB->isPseudo(Inst) || BC.MIB->isSystem(Inst)) {
            return false;
        }
    }

  return true;
}*/
bool PostLinkOutlining::isOutlinable(BinaryFunction &BF, 
                                     BinaryBasicBlock &BB, 
                                     const std::vector<MCInst> &Seq, 
                                     BinaryContext &BC) {
  // 1. 排除收缩包装 (Shrink-wrapping) 的区域
  // 通常 BOLT 会在 BF 或 BB 层面标记是否涉及到 shrink-wrapping

    if (!BF.isSimple()) {
    return false;
  }
  const auto *MIB = BC.MIB.get();
  const auto *MRI = BC.MRI.get();
const auto *MII = BC.MII.get();
unsigned i=0;
  for (const MCInst &Inst : Seq) {
    
    /*if (isInShrinkWrappedRegion(BB,Inst, BF)) {
    return false;
  }*/
  // 如果 i > 0（即不是序列的第一条指令），检查它是否是跳转目标
  if (i > 0) {
    if (BC.MIB->getTargetSymbol(Inst)) { // 检查指令是否带有目标符号
       return false;
    }
    // 更通用地，检查该指令在基本块中是否有 Label
    // 在 BOLT 中，通常通过 BC.getInstLabel(Inst) 获取关联符号
    if (MIB->getInstLabel(Inst)) { 
      return false;
    }
  }
    // 2. 排除跳转、分支、调用、返回
    if (MIB->isCall(Inst) || MIB->isBranch(Inst) || MIB->isReturn(Inst)) {
      return false;
    }

    // 3. 排除伪指令 (Pseudo) 和系统指令 (System/Privileged)
    if (MIB->isPseudo(Inst) ) {
      return false;
    }

    // 4. 排除对 PC, LR, FP 的使用 (读取、写入、计算)
    // 获取架构相关的特殊寄存器 ID
    
    MCPhysReg LR = BC.MRI->getRARegister(); // Link Register
    MCPhysReg FP = MIB->getFramePointer();  // Frame Pointer
    const MCInstrDesc &Desc = MII->get(Inst.getOpcode());
    for (unsigned i = 0; i < Inst.getNumOperands(); ++i) {
      const MCOperand &Op = Inst.getOperand(i);
      if (Op.isReg()) {
        MCPhysReg Reg = Op.getReg();
        
        // 检查是否是 LR 或 FP
        if (Reg == LR || Reg == FP) return false;
        
        // 检查寄存器别名（例如有些架构 FP 有多个名字）
        if (MRI->isSubRegisterEq(FP, Reg) || MRI->isSubRegisterEq(LR, Reg)) {
          return false;
        }
      }
      
      // 5. 排除涉及 PC 的隐式使用或计算 (如 PC-relative 寻址)
      // MIB 通常有专门的方法检测是否引用了 PC
      if (MIB->hasPCRelOperand(Inst)) {
        return false;
      }
    }

// 6. 隐式寄存器检查：排除隐式使用/定义 LR 和 FP
    // 在 LLVM 中，ImplicitUses/Defs 是指针，指向以 0 结尾的数组
    auto CheckImplicit = [&](ArrayRef<MCPhysReg> Regs) {
      for (MCPhysReg Reg : Regs) {
        if (MRI->isSubRegisterEq(LR, Reg) || (FP != 0 && MRI->isSubRegisterEq(FP, Reg))) {
          return true;
        }
      }
      return false;
    };

    if (CheckImplicit(Desc.implicit_uses()) || CheckImplicit(Desc.implicit_defs())) {
      return false;
    }
    i++;
  }

  return true;
}
void PostLinkOutlining::replaceSequencesWithCalls(
    BinaryContext &BC,
    std::vector<SequenceInstance> &Instances,
    BinaryFunction *OutlinedFunc){
  const auto *MIB = BC.MIB.get();
  
  for (int i = Instances.size() - 1; i >= 0; --i) {
    auto &seq = Instances[i];
    if (!seq.IsLabeled) continue;
    BinaryBasicBlock *BB = seq.BB;
    size_t SeqSize = seq.Insts.size();
    MCInst CallInst;
    BC.MIB->createCall(CallInst, OutlinedFunc->getSymbol(), BC.Ctx.get());
    auto CurrentIt = BB->begin() + seq.StartIndex;

    // 3. 删除旧的指令序列
    // eraseInstruction 会返回指向下一条指令的迭代器
    for (size_t i = 0; i < SeqSize; ++i) {
      CurrentIt = BB->eraseInstruction(CurrentIt);
    }

    // 4. 在删除的位置插入 Call 指令
    BB->insertInstruction(CurrentIt, CallInst);
    for (size_t j = i+1; j < Instances.size(); ++j) {
        // 跳过当前正在处理的这个序列
        if (i == (int)j) continue;
        
        auto &target = Instances[j];
        // 关键判断：同一个 BB，且目标序列在当前替换位置的【后面】
        if (target.BB == BB && target.StartIndex > seq.StartIndex) {
            target.StartIndex -= (SeqSize-1);
            
        }
        if(target.BB != BB)break;
    }
    seq.IsProcessed = true;
  }
  

}
// 3. 批量处理所有函数：遍历二进制上下文的函数集合
Error PostLinkOutlining::runOnFunctions(BinaryContext &BC) {
  unsigned index=0;
  errs()<<"begin "<<index<<"\n";
    for(int len = 32;len>=2;len--){
        auto seqs=getAllseqs(BC,len);
        uint64_t n=seqs.size();
        for(uint64_t i=0;i<n;i++){
            if(seqs[i].IsProcessed)continue;
            setLabel(seqs[i]);
            int frqcy=1;
            BinaryFunction *BF=seqs[i].BF;            // 所属函数
            BinaryBasicBlock *BB=seqs[i].BB;
            int b=seqs[i].StartIndex;
            for(uint64_t j=i+1;j<n;j++){
                if(seqs[j].IsProcessed)continue;
                if(areSequencesEquivalent(seqs[i].Insts,seqs[j].Insts,BC)&&!seqs[j].IsLabeled&&!(BF==seqs[j].BF&&BB==seqs[j].BB&&std::abs((int64_t)b - (int64_t)seqs[j].StartIndex) < (int64_t)len)){
                    setLabel(seqs[j]);
                    BF=seqs[j].BF;
                    BB=seqs[j].BB;
                    b=seqs[j].StartIndex;
                    frqcy++;
                }
            }
            BinaryFunction *OutlinedFunc;
            if(frqcy*len-frqcy-len-2>0){
              index++;
              errs()<<"outlined "<<index<<" len "<<len<<" frc "<<frqcy<<"\n";
              OutlinedFunc=createOutlinedFunction(BC,seqs[i].Insts,index);
              OutlinedFunc->setCodeSectionName(".text.injected");
              replaceSequencesWithCalls(BC,seqs,OutlinedFunc);
            }

            for(uint64_t j=0;j<n;j++){
              seqs[j].IsLabeled=false;
            }
        }
    }
  errs()<<"ended "<<index<<"\n";
  return Error::success();
}


}}