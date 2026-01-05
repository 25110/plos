#include "bolt/Core/BinaryFunctionCallGraph.h"
#include "bolt/Passes/BinaryPasses.h"
#include "bolt/Core/BinaryContext.h"
#include "bolt/Core/BinaryFunction.h"
#include "bolt/Passes/RegAnalysis.h"
namespace llvm {
namespace bolt {

class PostLinkOutlining : public BinaryFunctionPass {
public:
    struct SequenceInstance {
        std::vector<MCInst> Insts;      // 序列包含的指令指令流
        size_t Hash;                   // 预计算的哈希值，用于后续快速匹配
        BinaryFunction *BF;            // 所属函数
        BinaryBasicBlock *BB;          // 所属基本块
        BinaryBasicBlock::iterator It; // 序列在基本块中的起始位置迭代器
        bool IsLabeled = false;
        bool IsProcessed = false;
        // 用于判断重叠：该序列涵盖的指令地址范围（或在 BB 中的偏移）
        // 在 BOLT 中，我们可以简单地记录起始索引
        size_t StartIndex; 
    };
    // 1. 构造函数（接收“打印Pass日志”参数，遵循PPT示例）
    explicit PostLinkOutlining(const cl::opt<bool> &PrintPass)
        : BinaryFunctionPass(PrintPass) {}

    // 2. 唯一Pass名称（用于命令行调用和日志）
    const char *getName() const override { return "post-link-outlining"; }

    // 3. 核心执行方法：对整个二进制的函数批量处理（必重写）
    Error runOnFunctions(BinaryContext &BC) override;
    std::vector<SequenceInstance> getAllseqs(BinaryContext &BC, uint64_t len);
    void setLabel(SequenceInstance &Seq);
    bool isOutlinable(BinaryFunction &BF,BinaryBasicBlock &BB,const std::vector<MCInst> &Seq, BinaryContext &BC);
    bool hasMoreThanEightArgs(BinaryBasicBlock &BB, 
                                            BinaryBasicBlock::iterator CallIt, 
                                            BinaryContext &BC);
    bool isCalleeAccessingCallerStack(const MCInst &Inst,BinaryFunction *Callee, BinaryContext &BC);
    bool isInShrinkWrappedRegion(const BinaryBasicBlock &BB,const MCInst &Inst, const BinaryFunction &BF);
    bool PostLinkOutlining::areSequencesEquivalent(const std::vector<MCInst> &Seq1, 
                                               const std::vector<MCInst> &Seq2,
                                               const BinaryContext &BC);
    BinaryFunction* PostLinkOutlining::createOutlinedFunction(
    BinaryContext &BC, 
    const std::vector<MCInst> &Seq, 
    unsigned Index);
    void PostLinkOutlining::replaceSequencesWithCalls(
    BinaryContext &BC,
    std::vector<SequenceInstance> &Instances,
    BinaryFunction *OutlinedFunc);
private:

};

} }