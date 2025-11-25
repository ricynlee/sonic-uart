import tkinter as tk
from tkinter import scrolledtext, messagebox
import datetime
import random

class WeChatBubbleUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Sonic UART chat")
        self.root.geometry("400x600")
        self.root.configure(bg='#ECECEC')
        
        # 主框架
        self.main_frame = tk.Frame(root, bg='#ECECEC')
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # 聊天记录区域 - 使用Canvas实现滚动
        self.chat_canvas = tk.Canvas(
            self.main_frame, 
            bg='#F5F5F5',
            highlightthickness=0
        )
        self.chat_scrollbar = tk.Scrollbar(
            self.main_frame, 
            orient=tk.VERTICAL, 
            command=self.chat_canvas.yview
        )
        self.chat_canvas.configure(yscrollcommand=self.chat_scrollbar.set)
        
        self.chat_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.chat_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 聊天内容框架（放在Canvas上）
        self.chat_frame = tk.Frame(self.chat_canvas, bg='#F5F5F5')
        self.chat_window = self.chat_canvas.create_window(
            (0, 0), 
            window=self.chat_frame, 
            anchor=tk.NW,
            width=400
        )
        
        # 输入区域 - 使用Frame包装，确保按钮始终可见
        self.input_frame = tk.Frame(root, bg='#ECECEC')
        self.input_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # 使用Grid布局确保按钮在小窗口下也可见
        self.input_text = tk.Text(
            self.input_frame,
            height=3,
            font=('Microsoft YaHei', 10),
            relief=tk.FLAT,
            bg='white',
            wrap=tk.WORD
        )
        
        self.send_button = tk.Button(
            self.input_frame,
            text="发送",
            font=('Microsoft YaHei', 10),
            bg='#07C160',
            fg='white',
            relief=tk.FLAT,
            command=self.send_message,
            width=6
        )
        
        # 使用grid布局，输入框可伸缩，按钮固定宽度
        self.input_text.grid(row=0, column=0, sticky='ew', padx=(0, 5))
        self.send_button.grid(row=0, column=1, sticky='ew')
        
        # 设置列权重
        self.input_frame.grid_columnconfigure(0, weight=1)  # 输入框列可伸缩
        self.input_frame.grid_columnconfigure(1, weight=0)  # 按钮列固定宽度
        
        # 绑定事件
        self.root.bind('<Control-Return>', lambda event: self.send_message())
        self.input_text.bind('<Return>', self.handle_return_key)
        self.chat_frame.bind('<Configure>', self.on_frame_configure)
        self.chat_canvas.bind('<Configure>', self.on_canvas_configure)
        
        # 存储消息引用
        self.message_frames = []
        
        # 模拟一些初始消息
        self.init_messages()
    
    def on_frame_configure(self, event):
        """更新Canvas的滚动区域"""
        self.chat_canvas.configure(scrollregion=self.chat_canvas.bbox("all"))
    
    def on_canvas_configure(self, event):
        """调整内部frame的宽度以适应Canvas"""
        self.chat_canvas.itemconfig(self.chat_window, width=event.width)
    
    def handle_return_key(self, event):
        """处理回车键，Ctrl+Enter发送，普通回车换行"""
        if event.state & 0x4:  # Ctrl键被按下
            self.send_message()
            return "break"  # 阻止默认行为
        return None  # 允许正常换行
    
    def create_bubble(self, message, is_me=True):
        """创建聊天气泡"""
        # 主消息框架
        message_frame = tk.Frame(self.chat_frame, bg='#F5F5F5')
        self.message_frames.append(message_frame)
        
        # 时间显示框架
        time_frame = tk.Frame(message_frame, bg='#F5F5F5')
        time_str = datetime.datetime.now().strftime("%H:%M")
        time_label = tk.Label(
            time_frame,
            text=time_str,
            font=('Microsoft YaHei', 8),
            fg='#888888',
            bg='#F5F5F5'
        )
        time_label.pack(pady=2)
        
        # 内容框架（包含头像和气泡）
        content_frame = tk.Frame(message_frame, bg='#F5F5F5')
        
        # 气泡样式
        bubble_color = '#95EC69' if is_me else '#FFFFFF'
        bubble = tk.Label(
            content_frame,
            text=message,
            font=('Microsoft YaHei', 10),
            wraplength=250,
            justify=tk.LEFT,
            bg=bubble_color,
            fg='#000000',
            relief=tk.RAISED,
            bd=1,
            padx=12,
            pady=8
        )
        
        # 头像模拟
        avatar_color = '#FF6B6B' if is_me else '#4ECDC4'
        avatar_text = "此" if is_me else "彼"
        avatar = tk.Label(
            content_frame,
            text=avatar_text,
            font=('Microsoft YaHei', 8),
            width=3,
            height=1,
            bg=avatar_color,
            fg='white',
            relief=tk.FLAT
        )
        
        # 布局设置
        if is_me:
            # 我的消息：整体右对齐
            time_frame.pack(fill=tk.X, anchor=tk.E)
            content_frame.pack(fill=tk.X, anchor=tk.E)
            
            # 修正：头像在右，气泡也在右（气泡在头像左侧）
            avatar.pack(side=tk.RIGHT)  # 头像靠最右
            bubble.pack(side=tk.RIGHT, padx=(0, 5))  # 气泡在头像左边
            
        else:
            # 对方消息：整体左对齐
            time_frame.pack(fill=tk.X, anchor=tk.W)
            content_frame.pack(fill=tk.X, anchor=tk.W)
            
            # 对方消息：头像在左，气泡在右
            avatar.pack(side=tk.LEFT, padx=(0, 5))
            bubble.pack(side=tk.LEFT)
        
        return message_frame
    
    def add_message_to_chat(self, message, is_me=True):
        """添加消息到聊天区域"""
        message_frame = self.create_bubble(message, is_me)
        message_frame.pack(fill=tk.X, pady=5)
        
        # 滚动到底部
        self.chat_canvas.update_idletasks()
        self.chat_canvas.yview_moveto(1.0)
        self.chat_canvas.configure(scrollregion=self.chat_canvas.bbox("all"))
    
    def send_message(self):
        """发送消息"""
        message = self.input_text.get('1.0', tk.END).strip()
        if message:
            self.add_message_to_chat(message, is_me=True)
            self.input_text.delete('1.0', tk.END)
    
    def init_messages(self):
        """初始化一些示例消息"""
        messages = [
            ("对方发出的消息显示在此侧", False),
            ("我方发出的消息显示在此侧", True),
        ]
        
        for msg, is_me in messages:
            self.add_message_to_chat(msg, is_me)

def main():
    root = tk.Tk()
    app = WeChatBubbleUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
